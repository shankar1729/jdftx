/* Finite-difference tests for defect-only Coulomb electrostatics. GPL v3 or later. */
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/DefectCoulomb.h>
#include <electronic/SpeciesInfo.h>
#include <commands/parser.h>
#include <core/Util.h>
#include <cmath>

namespace
{
	void verify(double error, double tolerance, const char* label)
	{	logPrintf("%s: error %.6e (tolerance %.3e)\n", label, error, tolerance);
		if(!std::isfinite(error) || fabs(error)>tolerance) die("Test failed: %s\n", label);
	}
}

int main(int argc, char** argv)
{	InitParams ip("Check mixed Coulomb density/ionic derivatives against finite differences using a clean restart input.");
	initSystemCmdline(argc, argv, ip);
	{	Everything e;
		parse(readInputFile(ip.inputFilename), e);
		e.setup();
		auto correction = e.eVars.defectCoulomb;
		if(!correction) die("Test input must specify defect-coulomb and a clean restart.\n");
		e.eVars.n = e.eVars.calcDensity();
		ScalarField n0 = clone(e.eVars.get_nTot());
		ScalarFieldTilde source0 = clone(e.iInfo.getRhoIonBare());
		verify(correction->update(J(n0)), 1e-12, "Clean reference recovery");
		// Repeat evaluation to test that neither the frozen reference nor the ionic source is mutated.
		verify(correction->update(J(n0)), 1e-12, "Repeated clean recovery");
		verify(nrm2(e.iInfo.getRhoIonBare() - source0), 1e-12, "Ionic source preservation");

		// Two distinct, neutral density variations with strictly positive perturbed n.
		ScalarField p = clone(n0), q = clone(n0);
		const vector3<int>& S = e.gInfo.S;
		for(int ix=0; ix<S[0]; ix++) for(int iy=0; iy<S[1]; iy++) for(int iz=0; iz<S[2]; iz++)
		{	int i = (ix*S[1]+iy)*S[2]+iz;
			p->data()[i] *= sin(2*M_PI*ix/S[0]) + 0.3*sin(2*M_PI*iz/S[2]);
			q->data()[i] *= cos(2*M_PI*iy/S[1]) + 0.5*sin(2*M_PI*iz/S[2]);
		}
		p -= (sum(p)/sum(n0))*n0;
		q -= (sum(q)/sum(n0))*n0;
		ScalarField trial = n0 + 0.05*p;
		correction->update(J(trial));
		ScalarFieldTilde V = clone(correction->potential);
		double analytic = dot(J(q), O(V));
		const double h = 1e-4;
		double Eplus = correction->update(J(trial + h*q));
		double Eminus = correction->update(J(trial - h*q));
		verify((Eplus-Eminus)/(2*h) - analytic, 1e-10, "Density energy derivative");
		correction->update(J(trial));
		verify(nrm2(correction->potential - V), 1e-12, "Reference remains fixed after perturbations");

		// Verify actual insertion into the KS potential and energy, including JDFTx's dV normalization.
		e.eVars.n[0] = clone(trial);
		Energies corrected, baseline;
		e.eVars.EdensityAndVscloc(corrected);
		ScalarField correctedPotential = clone(e.eVars.Vscloc[0]);
		double Ecorr = correction->energy;
		e.eVars.defectCoulomb.reset();
		e.eVars.EdensityAndVscloc(baseline);
		verify(double(corrected.E) - double(baseline.E) - Ecorr, 1e-11, "Total energy insertion");
		verify(nrm2(correctedPotential - e.eVars.Vscloc[0] - e.gInfo.dV*I(V)), 1e-11, "KS potential insertion");
		e.eVars.defectCoulomb = correction;

		// Check explicit ionic derivative using the production pseudopotential force path.
		ScalarFieldTilde zero; initZero(zero, e.gInfo);
		ScalarFieldTilde nTrial = J(trial);
		correction->update(nTrial);
		auto sp = e.iInfo.species[0];
		ScalarFieldTilde coreZero = zero; // all optional gradients zero
		auto forces = sp->getLocalForces(zero, correction->ionicPotential, zero, coreZero, zero);
		vector3<> pos0 = sp->atpos[0];
		for(int dir=0; dir<3; dir++)
		{	sp->atpos[0] = pos0; sp->atpos[0][dir] += h; sp->sync_atpos();
			e.iInfo.update(e.ener);
			Eplus = correction->update(nTrial);
			sp->atpos[0] = pos0; sp->atpos[0][dir] -= h; sp->sync_atpos();
			e.iInfo.update(e.ener);
			Eminus = correction->update(nTrial);
			verify((Eplus-Eminus)/(2*h) + forces[0][dir], 2e-7, "Ionic energy derivative (lattice coordinates)");
		}
		sp->atpos[0] = pos0; sp->sync_atpos();
		e.iInfo.update(e.ener);
		if(e.coulombParams.embed)
		{	// An independent check against JDFTx's public point-charge operators:
			// D*(deltaN+S*deltaIon) must equal D*deltaN + (Kiso^R-Kslab^R)*deltaIon.
			sp->atpos[0][0] += 0.01; sp->sync_atpos(); e.iInfo.update(e.ener);
			correction->update(nTrial);
			CoulombParams iso = e.coulombParams;
			iso.geometry = CoulombParams::Isolated;
			iso.exchangeRegularization = CoulombParams::None; iso.omegaSet.clear();
			if(correction->hasCenter) iso.embedCenter = correction->center;
			auto kernel = iso.createCoulomb(e.gInfo, " for independent point-charge check");
			ScalarFieldTilde dn = J(trial-n0), di = e.iInfo.getRhoIonBare()-source0;
			ScalarFieldTilde expected = J(I((*kernel)(dn) - (*e.coulomb)(dn)
				+ (*kernel)(di, Coulomb::PointChargeRight) - (*e.coulomb)(di, Coulomb::PointChargeRight)));
			verify(nrm2(correction->potential-expected), 1e-9, "Embedded PointChargeRight equivalence");
			sp->atpos[0] = pos0; sp->sync_atpos();
		}
		logPrintf("All defect Coulomb derivative tests passed.\n");
	}
	finalizeSystem();
	return 0;
}
