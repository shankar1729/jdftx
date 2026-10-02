/* Defect-only mixed-boundary electrostatics. Licensed under GPL v3 or later. */
#include <electronic/DefectCoulomb.h>
#include <electronic/Everything.h>
#include <electronic/SpeciesInfo.h>
#include <core/ScalarFieldIO.h>
#include <cstdint>
#include <cmath>
#include <cstring>

namespace
{
	// Fingerprint the actual pseudopotential bytes, independent of its location.
	unsigned long long fingerprint(const string& filename)
	{	FILE* fp = fopen(filename.c_str(), "rb");
		if(!fp) die("Cannot read pseudopotential '%s' for defect reference validation.\n", filename.c_str());
		uint64_t hash = UINT64_C(14695981039346656037);
		unsigned char buffer[8192];
		size_t count;
		while((count = fread(buffer, 1, sizeof(buffer), fp)))
			for(size_t i=0; i<count; i++) { hash ^= buffer[i]; hash *= UINT64_C(1099511628211); }
		if(ferror(fp)) die("Error reading pseudopotential '%s'.\n", filename.c_str());
		fclose(fp);
		return hash;
	}

	void check(bool valid, const char* what)
	{	if(!valid) die("Invalid or incompatible defect Coulomb reference: %s.\n", what);
	}

	bool same(double x, double y)
	{	return std::isfinite(x) && std::isfinite(y) && fabs(x-y) <= 1e-12*std::max(1., fabs(y));
	}

	void checkReferenceMode(const Everything& e)
	{	check(e.coulombParams.geometry == CoulombParams::Slab, "requires coulomb-interaction Slab");
		check(!e.coulombParams.embedFluidMode, "fluid-mode Coulomb embedding is not supported");
		check(e.eVars.fluidParams.fluidType == FluidNone, "requires a vacuum calculation");
		check(!e.iInfo.ljOverride, "requires electronic structure, not pair-potential override");
		check(!e.eVars.rhoExternal && !e.coulombParams.Efield.length_squared(), "external charges/fields are not supported");
	}
}

void DefectCoulomb::saveReference(const Everything& e, const string& filename)
{	checkReferenceMode(e);
	check(!e.eVars.defectCoulomb, "generate the clean reference without defect-coulomb");
	check(bool(e.eVars.n[0]), "electronic density has not been initialized");
	ScalarField n = e.eVars.get_nTot();
	const ScalarFieldTilde& ion = e.iInfo.getRhoIonBare();
	double q = e.gInfo.detR * (J(n)->getGzero() + ion->getGzero());
	check(std::isfinite(q) && fabs(q) <= 1e-5, "clean reference must be neutral");
	if(!mpiWorld->isHead()) return;
	FILE* fp = fopen(filename.c_str(), "wb");
	if(!fp) die("Cannot write defect Coulomb reference '%s'.\n", filename.c_str());
	fprintf(fp, "JDFTX_DEFECT_REFERENCE 2\n");
	fprintf(fp, "%d %d %d\n", e.gInfo.S[0], e.gInfo.S[1], e.gInfo.S[2]);
	for(int i=0; i<3; i++) for(int j=0; j<3; j++) fprintf(fp, "%.17g\n", e.gInfo.R(i,j));
	fprintf(fp, "%.17g %.17g %d %d %d %.17g\n", e.cntrl.Ecut, e.cntrl.EcutRho,
		e.coulombParams.iDir, int(e.eInfo.spinType), int(e.eInfo.smearingType), e.eInfo.smearingWidth);
	int nSpecies = 0;
	for(auto sp: e.iInfo.species) if(sp->atpos.size()) nSpecies++;
	fprintf(fp, "%d\n", nSpecies);
	for(auto sp: e.iInfo.species) if(sp->atpos.size())
		fprintf(fp, "%s %.17g %016llx\n", sp->name.c_str(), sp->Z, fingerprint(sp->potfilename));
	fprintf(fp, "%d %.17g %.17g %.17g %.17g\n", int(e.coulombParams.embed),
		e.coulombParams.embedCenter[0], e.coulombParams.embedCenter[1], e.coulombParams.embedCenter[2], e.coulombParams.ionMargin);
	fprintf(fp, "DATA\n");
	// Preserve all point-ion Fourier coefficients: embedded range separation must
	// precede the real-grid projection, particularly near Nyquist frequencies.
	saveRawBinary(n, fp);
	saveRawBinary(ion, fp);
	if(fclose(fp)) die("Failed to close defect Coulomb reference '%s'.\n", filename.c_str());
}

void DefectCoulomb::setup(const Everything& everything, bool fromFixedDensity)
{	e = &everything;
	checkReferenceMode(*e);
	if(fromFixedDensity)
		check(e->cntrl.fixed_H && e->eVars.nFilenamePattern.length(), "one-shot mode requires a supplied electron density");
	else
		check(e->cntrl.scf && !e->cntrl.fixed_H, "requires electronic-scf with a variable density");
	check(e->symm.mode == SymmetriesNone, "use symmetries none to preserve the fixed reference");
	check(!e->ionicMinParams.nIterations && !e->latticeMinParams.nIterations && !e->ionicDynParams.nSteps,
		"this implementation requires fixed ions and lattice");
	check(!e->iInfo.computeStress && !e->vibrations && !e->pertInfo.solverParams.nIterations,
		"stress, vibrations and perturbation calculations are not supported");
	check(std::isnan(e->eInfo.mu), "fixed chemical potential is not supported; use neutral fixed electron number");
	check(fabs(e->eInfo.nElectrons - e->iInfo.getZtot()) <= chargeTolerance, "defect system must be neutral");
	logPrintf("\nReading fixed clean reference '%s' ... ", referenceFilename.c_str()); logFlush();
	FILE* fp = fopen(referenceFilename.c_str(), "rb");
	if(!fp) die("Cannot read defect Coulomb reference '%s'.\n", referenceFilename.c_str());
	char magic[64]; int version;
	check(fscanf(fp, "%63s %d", magic, &version)==2 && !strcmp(magic, "JDFTX_DEFECT_REFERENCE") && (version==1 || version==2), "file format/version");
	vector3<int> S;
	check(fscanf(fp, "%d %d %d", &S[0], &S[1], &S[2])==3 && S==e->gInfo.S, "FFT grid mismatch");
	for(int i=0; i<3; i++) for(int j=0; j<3; j++)
	{	double value;
		check(fscanf(fp, "%lf", &value)==1 && same(value, e->gInfo.R(i,j)), "lattice mismatch");
	}
	double Ecut, EcutRho, width; int iDir, spin, smearing;
	check(fscanf(fp, "%lf %lf %d %d %d %lf", &Ecut, &EcutRho, &iDir, &spin, &smearing, &width)==6, "settings header");
	check(same(Ecut, e->cntrl.Ecut) && same(EcutRho, e->cntrl.EcutRho), "plane-wave cutoff mismatch");
	check(iDir==e->coulombParams.iDir && spin==int(e->eInfo.spinType), "slab direction/spin mismatch");
	check(smearing==int(e->eInfo.smearingType) && same(width, e->eInfo.smearingWidth), "smearing mismatch");
	int nSpecies;
	check(fscanf(fp, "%d", &nSpecies)==1 && nSpecies>=0 && nSpecies<=int(e->iInfo.species.size()), "host species count");
	std::set<string> seen;
	for(int i=0; i<nSpecies; i++)
	{	char name[256]; double Z; unsigned long long hash;
		check(fscanf(fp, "%255s %lf %llx", name, &Z, &hash)==3, "host species header");
		check(seen.insert(name).second, "duplicate host species");
		bool found = false;
		for(auto sp: e->iInfo.species) if(sp->name==name)
		{	check(same(Z, sp->Z) && hash==fingerprint(sp->potfilename), "host pseudopotential mismatch");
			found = true;
		}
		check(found, "host species missing from defect system");
	}
	if(version==2)
	{	int embed; vector3<> referenceCenter; double margin;
		check(fscanf(fp, "%d %lf %lf %lf %lf", &embed, &referenceCenter[0], &referenceCenter[1], &referenceCenter[2], &margin)==5,
			"embedding header");
		check(embed==int(e->coulombParams.embed), "embedding mode mismatch");
		if(embed)
		{	for(int k=0; k<3; k++)
			{	double offset = referenceCenter[k] - e->coulombParams.embedCenter[k];
				check(std::isfinite(offset) && fabs(offset-round(offset))<=1e-12, "embedding center mismatch");
			}
			check(same(margin, e->coulombParams.ionMargin), "embedding ionic margin mismatch");
		}
	}
	else check(!e->coulombParams.embed, "version-1 references do not record embedding; regenerate the clean reference");
	check(fscanf(fp, "%63s", magic)==1 && !strcmp(magic, "DATA") && fgetc(fp)=='\n', "data marker");
	long dataOffset = ftell(fp);
	check(dataOffset>=0 && fileSize(referenceFilename.c_str()) == dataOffset
		+ off_t(e->gInfo.nr)*off_t(sizeof(double))
		+ (version==1 ? off_t(e->gInfo.nr)*off_t(sizeof(double)) : off_t(e->gInfo.nG)*off_t(sizeof(complex))), "density payload length");
	ScalarField n = ScalarFieldData::alloc(e->gInfo);
	loadRawBinary(n, fp);
	ionReference = ScalarFieldTildeData::alloc(e->gInfo);
	if(version==1)
	{	ScalarField ion = ScalarFieldData::alloc(e->gInfo);
		loadRawBinary(ion, fp); ionReference = J(ion);
	}
	else loadRawBinary(ionReference, fp);
	fclose(fp);
	for(int i=0; i<e->gInfo.nr; i++) check(std::isfinite(n->data()[i]), "non-finite electronic reference density");
	for(int i=0; i<e->gInfo.nG; i++)
		check(std::isfinite(ionReference->data()[i].real()) && std::isfinite(ionReference->data()[i].imag()), "non-finite ionic reference density");
	nReference = J(n);
	if(!e->coulombParams.embed) ionReference = J(I(ionReference));
	double q = e->gInfo.detR * (nReference->getGzero() + ionReference->getGzero());
	check(std::isfinite(q) && fabs(q)<=chargeTolerance, "non-neutral clean reference");
	logPrintf("done (net reference charge %+.3e e)\n", q);
	isolatedParams = e->coulombParams;
	isolatedParams.geometry = CoulombParams::Isolated;
	if(hasCenter)
	{	check(e->coulombParams.embed, "defect-coulomb-center requires coulomb-truncation-embed");
		for(int k=0; k<3; k++) check(std::isfinite(center[k]), "non-finite defect center");
		isolatedParams.embedCenter = center;
	}
	isolatedParams.exchangeRegularization = CoulombParams::None;
	isolatedParams.omegaSet.clear();
	isolatedParams.computeStress = false;
	isolated = isolatedParams.createCoulomb(e->gInfo, " for localized defect charge only");
	check(same(isolated->getIonWidth(), e->coulomb->getIonWidth()), "inconsistent ionic range separation");
	if(e->coulombParams.embed)
		logPrintf("Defect Coulomb: embedded slab and isolated operators; ionic Gaussian width %.12g bohr.\n", isolated->getIonWidth());
	logPrintf("Defect Coulomb: fixed reference; isolated WS correction to slab electrostatics.\n");
}

double DefectCoulomb::update(const ScalarFieldTilde& nTilde)
{	ScalarFieldTilde deltaN = nTilde - nReference;
	ScalarFieldTilde deltaIon = e->coulombParams.embed ? e->iInfo.getRhoIonBare() - ionReference
		: J(I(e->iInfo.getRhoIonBare())) - ionReference;
	deltaRho = deltaN + J(I(deltaIon));
	charge = e->gInfo.detR * deltaRho->getGzero();
	if(!std::isfinite(charge) || fabs(charge)>chargeTolerance)
		die("Defect Coulomb density difference is not neutral: %+.12g e (tolerance %.3g).\n", charge, chargeTolerance);
	// Do not discard G=0: the finite truncated-kernel value fixes the potential gauge.
	// Project the output too: P(Kiso-Kslab)P is self-adjoint on the real FFT grid.
	// Both embedded operators use the same ionic short-range kernel, which cancels
	// in their difference. Their PointChargeRight difference is therefore D*S,
	// where D=Kiso-Kslab on smooth sources and S is the ionic Gaussian filter.
	// The variational energy is 1/2 (deltaN+S*deltaIon) D (deltaN+S*deltaIon).
	// Its electron gradient is V=D*(deltaN+S*deltaIon), and ion gradient is S*V.
	double width = isolated->getIonWidth();
	effectiveSource = deltaN + J(I(gaussConvolve(deltaIon, width)));
	potential = J(I((*isolated)(effectiveSource) - (*e->coulomb)(effectiveSource)));
	ionicPotential = gaussConvolve(potential, width);
	energy = 0.5 * dot(effectiveSource, O(potential));
	return energy;
}

void DefectCoulomb::report() const
{	logPrintf("DefectCoulomb: Qdefect: %+.6e e  Ecorr: %+.12e Eh\n", charge, energy);
}
