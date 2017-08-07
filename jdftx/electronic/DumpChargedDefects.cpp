#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/Dump_internal.h>
#include <core/WignerSeitz.h>
#include <core/Operators.h>
#include <core/ScalarFieldIO.h>
#include <core/LatticeUtils.h>
#include <core/Coulomb_internal.h>
#include <fluid/PCM.h>
#include <gsl/gsl_sf.h>

//-------------------------- Slab epsilon ----------------------------------

inline void planarAvg_sub(size_t iStart, size_t iStop, const vector3<int>& S, int iDir, complex* data)
{	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	THREAD_halfGspaceLoop( if(iG[jDir] || iG[kDir]) data[i] = 0.; )
}
inline void planarAvg(ScalarFieldTilde& X, int iDir)
{	threadLaunch(planarAvg_sub, X->gInfo.nG, X->gInfo.S, iDir, X->data());
}
inline ScalarField getPlanarAvg(const ScalarField& X, int iDir)
{	ScalarFieldTilde Xtilde = J(X);
	planarAvg(Xtilde, iDir);
	return I(Xtilde);
}
inline void planarAvg(ScalarField& X, int iDir)
{	X = getPlanarAvg(X, iDir);
}

//Multiply by phase factor, that when averaged extracts the prefactor of -sin(2 pi x)/(2pi) (the form used in Coulomb::getEfieldPotential) along direction jDir
inline void sinMultiply_sub(size_t iStart, size_t iStop, const vector3<int>& S, int jDir, const vector3<>& xCenter, double* data)
{	double invSj = 1./S[jDir];
	THREAD_rLoop(
		double xj = invSj * iv[jDir] - xCenter[jDir];
		data[i] *= (-4*M_PI) * sin(2*M_PI*xj);
	)
}
inline void sinMultiply(ScalarField& X, int jDir, const vector3<>& xCenter)
{	threadLaunch(sinMultiply_sub, X->gInfo.nr, X->gInfo.S, jDir, xCenter, X->data());
}

inline void fixBoundary_sub(size_t iStart, size_t iStop, const vector3<int>& S, int iDir, int iBoundary, double* eps)
{	iBoundary = positiveRemainder(iBoundary, S[iDir]);
	THREAD_rLoop(
		int dist = iv[iDir] - iBoundary;
		if(dist*2 > S[iDir]) dist -= S[iDir];
		if(dist*2 < -S[iDir]) dist += S[iDir];
		if(abs(dist)<=2) eps[i] = 1.;
	)
}
inline void fixBoundarySmooth(ScalarField& epsInv, int iDir, const vector3<>& xCenter, double sigma)
{	const GridInfo& gInfo = epsInv->gInfo;
	//Fix values of epsilon near truncation boundary to 1:
	int iBoundary = int(round((xCenter[iDir]+0.5) * gInfo.S[iDir]));
	threadLaunch(fixBoundary_sub, gInfo.nr, gInfo.S, iDir, iBoundary, epsInv->data());
	//Smooth:
	epsInv = I(gaussConvolve(J(epsInv), sigma));
}

void SlabEpsilon::dump(const Everything& e, ScalarField d_tot) const
{	string fname = e.dump.getFilename("slabEpsilon");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	
	//Find change in d_tot and applied fields:
	ScalarField dDiff;
	{	ScalarField d_totRef(ScalarFieldData::alloc(e.gInfo));
		loadRawBinary(d_totRef, dtotFname.c_str());
		dDiff = d_tot - d_totRef;
	}
	vector3<> Ediff = e.coulombParams.Efield - Efield;
	assert(Ediff.length_squared());
	
	//Calculate inverse of epsilon along truncated direction:
	int iDir = e.coulombParams.iDir;
	double h = e.gInfo.h[iDir].length();
	vector3<> zHat = e.gInfo.R.column(iDir);
	zHat *= 1./zHat.length(); //unit vector along slab normal
	double EzDiff = dot(zHat, Ediff);
	ScalarField epsInv_z;
	if(EzDiff)
	{	//Calculate field using central-difference derivative:
		ScalarFieldTilde tPlus(ScalarFieldTildeData::alloc(e.gInfo)), tMinus(ScalarFieldTildeData::alloc(e.gInfo));
		initTranslation(tPlus, e.gInfo.h[iDir]);
		initTranslation(tMinus, -e.gInfo.h[iDir]);
		ScalarFieldTilde epsInvTilde = (1./(EzDiff * 2.*h)) * (tPlus - tMinus) * J(dDiff);
		planarAvg(epsInvTilde, iDir);
		epsInv_z = I(epsInvTilde);
		fixBoundarySmooth(epsInv_z, iDir, e.coulomb->xCenter, sigma);
	}
	
	//Write file:
	FILE* fp = fopen(fname.c_str(), "w");
	fprintf(fp, "#distance[bohr]  epsilon_normal  epsilon_||\n");
	vector3<int> iR;
	double* epsInv_zData = epsInv_z ? epsInv_z->data() : 0;
	for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
	{	double z = h*iR[iDir];
		size_t i = e.gInfo.fullRindex(iR);
		fprintf(fp, "%lf %lf\n", z,
			(epsInv_z ? 1./epsInv_zData[i] : NAN) );
	}
	fclose(fp);
	logPrintf("done\n"); logFlush();
}

//----------- Bulk epsilon ------

void BulkEpsilon::dump(const Everything& e, ScalarField d_tot) const
{	//Find change in d_tot and applied fields:
	ScalarField dDiff;
	{	ScalarField d_totRef(ScalarFieldData::alloc(e.gInfo));
		loadRawBinary(d_totRef, dtotFname.c_str());
		dDiff = d_tot - d_totRef;
	}
	vector3<> Ediff = e.coulombParams.Efield - Efield;
	assert(Ediff.length_squared());

	//Calculate inverse of epsilon along parallel directions:
	vector3<> RT_Ediff = e.gInfo.RT * Ediff;
	double wSum = 0., epsInv = 0.;
	for(int jDir=0; jDir<3; jDir++)
	{	double w = fabs(RT_Ediff[jDir])/(e.gInfo.R.column(jDir).length() * Ediff.length());
		if(w < symmThreshold) continue;
		ScalarField epsInvCur = (1./RT_Ediff[jDir]) * dDiff;
		sinMultiply(epsInvCur, jDir, e.coulombParams.embedCenter); //this extracts the sin component upon unit-cell averaging below
		epsInv += w * integral(epsInvCur)/e.gInfo.detR; //average over unit cell
		wSum += w;
	}
	epsInv *= (1./wSum); //Now contains the direction-averaged response (for whichever components are available)
	logPrintf("bulkEpsilon = %lg\n", 1./epsInv);
}

//-------------------------- Charged defects ----------------------------------

struct SlabPeriodicSolver : public LinearSolvable<ScalarFieldTilde>
{	int iDir; //truncated direction
	VectorField epsMinus1; //anisotropic epsilon - 1
	const ScalarField& kappaSq;
	const GridInfo& gInfo;
	RealKernel Ksqrt, Kinv;
	const ScalarField epsInv;
	double K0; //G=0 component of slab kernel
	
	static inline void setKernels_sub(size_t iStart, size_t iStop, const GridInfo* gInfo, int iDir, double kRMS, double* Ksqrt, double* Kinv, bool periodicCoulomb)
	{	const vector3<int>& S = gInfo->S;
		CoulombSlab_calc slabCalc(iDir, 0.5*gInfo->R.column(iDir).length());
		THREAD_halfGspaceLoop
		(	if(periodicCoulomb)
				Kinv[i] = gInfo->GGT.metric_length_squared(iG); //no truncation; image separation is due to fluid response
			else
			{	double K = slabCalc(iG, gInfo->GGT) / (4*M_PI);
				Kinv[i] = (fabs(K)>1e-12) ? 1./K : 0.;
			}
			Ksqrt[i] = (Kinv[i] || kRMS) ? 1./sqrt(fabs(Kinv[i]) + kRMS*kRMS) : 0.;
		)
	}
	
	SlabPeriodicSolver(int iDir, const ScalarField& epsPerp, const ScalarField& epsPar, const ScalarField& kappaSq, bool periodicCoulomb)
	: iDir(iDir), kappaSq(kappaSq), gInfo(epsPerp->gInfo), Ksqrt(gInfo), Kinv(gInfo), epsInv(3.*inv(epsPerp+2.*epsPar))
	{
		for(int jDir=0; jDir<3; jDir++)
			epsMinus1[jDir] = (jDir==iDir) ? (epsPerp-1.) : (epsPar-1.);
		
		double kRMS = sqrt(integral(kappaSq)*3./integral(epsPerp+2.*epsPar)); //average Debye length of unit cell (for preconditioner)
		threadLaunch(setKernels_sub, gInfo.nG, &gInfo, iDir, kRMS, Ksqrt.data(), Kinv.data(), periodicCoulomb);
		nullToZero(state, gInfo);
	}
	
	ScalarFieldTilde hessian(const ScalarFieldTilde& phiTilde) const
	{	ScalarFieldTilde rhoTilde = -(Kinv * phiTilde); //vacuum term
		rhoTilde += divergence(J(epsMinus1 * I(gradient(phiTilde))));  //dielectric term
		rhoTilde -= J(kappaSq * I(phiTilde)); //Debye screening term
		return (-1./(4*M_PI)) * rhoTilde;
	}
	
	ScalarFieldTilde precondition(const ScalarFieldTilde& rTilde) const
	{	return Ksqrt*(J(epsInv*I(Ksqrt*rTilde)));
	}

	double getEnergy(const ScalarFieldTilde& rho, ScalarFieldTilde& phi)
	{	MinimizeParams mp;
		mp.nDim = gInfo.nr;
		mp.knormThreshold = 1e-11;
		mp.fpLog = globalLog;
		mp.linePrefix = "\tSlabPeriodicCG: ";
		mp.energyFormat = "%+.15lf";

		zeroNyquist((ScalarFieldTilde&)rho);
		solve(rho, mp);
		
		phi = state;
		return 0.5*dot(phi,O(rho));
	}
};


struct CylindricalPoisson
{
	int NZ; //number of grid points along truncated direction
	double L; //length along truncated direction
	double epsilonBulk, kappaSqBulk; //bulk response
	matrix epsParTilde, GGepsKappaSq; diagMatrix G; //Fourier space operators
	
	CylindricalPoisson(int iDir, const ScalarField& epsPerpSlab, const ScalarField& epsParSlab, const ScalarField& kappaSqSlab)
	{
		//Extract grid dimensions and bulk response:
		const GridInfo& gInfo = epsPerpSlab->gInfo;
		NZ = gInfo.S[iDir];
		L = gInfo.R.column(iDir).length();
		vector3<int> iRbulk; iRbulk[iDir] = NZ/2;
		size_t iBulk = gInfo.fullRindex(iRbulk);
		epsilonBulk = (epsPerpSlab->data()[iBulk] + 2.*epsParSlab->data()[iBulk])/3; //should be isotropic anyway
		kappaSqBulk = kappaSqSlab->data()[iBulk];
		
		//Initialize Fourier space operators along truncated direction:
		G.resize(NZ);
		for(int iZ=0; iZ<NZ; iZ++)
			G[iZ] = (2*M_PI/L) * (iZ<NZ/2 ? iZ : iZ-NZ);
		complexScalarFieldTilde epsPerpSlabTilde = J(Complex(epsPerpSlab - epsilonBulk));
		complexScalarFieldTilde epsParSlabTilde  = J(Complex(epsParSlab  - epsilonBulk ));
		complexScalarFieldTilde kappaSqSlabTilde = J(Complex(kappaSqSlab - kappaSqBulk));
		std::vector<complex> epsPerpDiagTilde(NZ), epsParDiagTilde(NZ), kappaSqDiagTilde(NZ);
		for(int iZ=0; iZ<NZ; iZ++)
		{	vector3<int> iG; iG[iDir]=iZ;
			size_t iSlab = gInfo.fullRindex(iG);
			epsPerpDiagTilde[iZ] = epsPerpSlabTilde->data()[iSlab];
			epsParDiagTilde[iZ]  = epsParSlabTilde->data()[iSlab];
			kappaSqDiagTilde[iZ] = kappaSqSlabTilde->data()[iSlab];
		}
		epsParTilde.init(NZ,NZ);
		GGepsKappaSq.init(NZ,NZ);
		for(int iZ=0; iZ<NZ; iZ++)
			for(int jZ=0; jZ<NZ; jZ++)
			{	int kZ = (iZ - jZ);
				if(kZ<0) kZ += NZ; //wrap kZ to [0,NZ)
				epsParTilde.set(iZ,jZ, epsParDiagTilde[kZ]);
				GGepsKappaSq.set(iZ,jZ, G[iZ]*epsPerpDiagTilde[kZ]*G[jZ] + kappaSqDiagTilde[kZ]);
			}
	}
    
	//Integrand
	double integrand(double k, double sigma, const matrix& OgTilde) const
	{	matrix KinvTot = zeroes(NZ,NZ);
		//Set truncated Greens function
		double alpha = sqrt(k*k + kappaSqBulk/epsilonBulk);
		double expMhlfAlphaL = exp(-0.5*alpha*L), cosMhlfGL = 1.;
		for(int iZ=0; iZ<NZ; iZ++)
		{	KinvTot.set(iZ,iZ, L*epsilonBulk*(alpha*alpha + G[iZ]*G[iZ])/(1. - expMhlfAlphaL*cosMhlfGL));
			cosMhlfGL = -cosMhlfGL; //since Gn L = 2 n pi
		}
		//Add inhomogeneous screening terms:
		KinvTot += L * (GGepsKappaSq + (k*k)*epsParTilde);
		//Calculate self-energy at k:
		double Uk = trace(dagger(OgTilde) * invApply(KinvTot, OgTilde)).real();
		return k * exp(-std::pow(k*sigma,2)) * Uk;
	}
	struct IntegrandParams { double sigma; const matrix* OgTilde; const CylindricalPoisson* cp; };
	static double integrand_wrapper(double k, void* params) //wrapper for GSL integration routine
	{	const IntegrandParams& ip = *((const IntegrandParams*)params);
		return ip.cp->integrand(k, ip.sigma, *(ip.OgTilde));
	}
	
	//Calculate self energy of Gaussian with norm q and width sigma centered at z0
	double getEnergy(double q, double sigma, double z0) const
	{	//Source term:
		matrix OgTilde(NZ,1);
		for(int iZ=0; iZ<NZ; iZ++)
			OgTilde.set(iZ,0, exp(-0.5*std::pow(G[iZ]*sigma,2))*cis(-G[iZ]*z0)); //note O cancels 1./L in derivation
		size_t wsSize = 1024;
		gsl_integration_workspace* ws = gsl_integration_workspace_alloc(wsSize);
		IntegrandParams ip = { sigma, &OgTilde, this };
		gsl_function f;
		f.function = integrand_wrapper;
		f.params = &ip;
		double integral, intErr;
		gsl_integration_qagiu(&f, 0., 1e-12, 1e-12, wsSize, ws, &integral, &intErr); //Calculate \int_0^infty dk integrand(k)
		gsl_integration_workspace_free(ws);
		return (q*q) * integral;
	}
};


void ChargedDefect::dump(const Everything& e, ScalarField d_tot) const
{	logPrintf("Calculating charged defect correction:\n"); logFlush();
	if(!center.size())
		die("\tNo model charges specified (using command charged-defect).\n");
	
	//Construct model charge on plane-wave grid:
	ScalarFieldTilde rhoModel;
	double qTot = 0.;
	for(const Center& cdc: center)
	{	ScalarFieldTilde trans(ScalarFieldTildeData::alloc(e.gInfo));
		initTranslation(trans, e.gInfo.R * cdc.pos);
		rhoModel += gaussConvolve((cdc.q/e.gInfo.detR)*trans, cdc.sigma);
		qTot += cdc.q;
	}
	
	//Find center of defects for alignment calculations:
	logSuspend(); WignerSeitz ws(e.gInfo.R); logResume();
	vector3<> pos0 = geometry==CoulombParams::Periodic ? center[0].pos : e.coulomb->xCenter;
	vector3<> posMean = pos0;
	for(const Center& cdc: center)
		posMean += (1./center.size()) * ws.restrict(cdc.pos - pos0);
	
	//Read reference Dtot and calculate electrostatic potential difference within DFT:
	ScalarField d_totRef(ScalarFieldData::alloc(e.gInfo));
	loadRawBinary(d_totRef, dtotFname.c_str());
	ScalarField Vdft = d_tot - d_totRef; //electrostatic potential of defect from DFT
	
	//Calculate isolated and periodic self-energy (and potential) of model charge
	ScalarField Vmodel; double Emodel=0., EmodelIsolated=0.;
	switch(geometry)
	{	case CoulombParams::Periodic: //Bulk defect
		{	//Periodic potential and energy:
			ScalarFieldTilde dModel = (*e.coulomb)(rhoModel) * (1./bulkEps); //assuming uniform dielectric
			Emodel = 0.5*dot(rhoModel, O(dModel));
			Vmodel = I(dModel);
			//Isolated energy:
			for(const Center& cdc: center)
				EmodelIsolated += 0.5*std::pow(cdc.q,2)/(cdc.sigma*sqrt(M_PI)*bulkEps); //self energy of Gaussian
			break;
		}
		case CoulombParams::Slab: //Surface defect
		{	bool truncated = (e.coulombParams.geometry==CoulombParams::Slab);
			if(truncated && !e.coulombParams.embed)
				die("\tCoulomb truncation must be embedded for charged-defect correction in slab geometry.\n");
			if(truncated) rhoModel = e.coulomb->embedExpand(rhoModel); //switch to embedding grid
			
			//Create dielectric model for slab:
			ScalarFieldArray epsSlab; nullToZero(epsSlab, e.gInfo, 2); //2 components = perp, par
			double* epsSlabData[2] = { epsSlab[0]->data(), epsSlab[1]->data() };
			std::ifstream ifs(slabEpsFname.c_str());
			if(!ifs.is_open()) die("\tCould not open slab dielectric model file '%s' for reading.\n", slabEpsFname.c_str());
			string commentLine; getline(ifs, commentLine); //get and ignore comment line
			vector3<int> iR;
			double h = e.gInfo.h[iDir].length();
			for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
			{	//Read data from line:
				string epsLine; getline(ifs, epsLine);
				istringstream iss(epsLine);
				double d;
				double &epsPerp = epsSlabData[0][e.gInfo.fullRindex(iR)];
				double &epsPar = epsSlabData[1][e.gInfo.fullRindex(iR)];
				iss >> d >> epsPerp;
				iss >> epsPar;
				if(iss.fail()) epsPar = epsPerp; //seperate parallel value not available
				//Check grid point:
				double dExpected = h*iR[iDir];
				if(fabs(d-dExpected) > symmThreshold)
					die("\tGeometry mismatch in '%s': expecting distance %lg on line %d; found %lg instead.\n",
						slabEpsFname.c_str(), dExpected, iR[iDir]+2, d);
			}
			for(int dir=0; dir<2; dir++)
			{	ScalarFieldTilde epsSlabMinus1tilde = J(epsSlab[dir]) * (e.gInfo.nr / e.gInfo.S[iDir]); //multiply by number of points per plane (to account for average below)
				epsSlabMinus1tilde->setGzero(epsSlabMinus1tilde->getGzero() - 1.); //subtract 1
				planarAvg(epsSlabMinus1tilde, iDir); //now contains a planarly-uniform version of epsSlab-1
				if(truncated)
					epsSlabMinus1tilde = e.coulomb->embedExpand(epsSlabMinus1tilde); //switch to embedding grid (note embedding eps-1 because it is zero far away)
				epsSlab[dir] = 1. + I(epsSlabMinus1tilde);
			}
			
			//Include solvation model dielectric / screening, if present:
			ScalarField kappaSqSlab; nullToZero(kappaSqSlab, epsSlab[0]->gInfo);
			double epsilonBulk = 1., kappaSqBulk = 0.;
			if(e.eVars.fluidSolver)
			{	if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
					logPrintf("WARNING: charged-defect-correction does not support ClassicalDFT; ignoring fluid response\n");
				else
				{	//Fluid is a PCM: use cavity shape to update epsilon, kappaSq
					//(approx. fluid as linear and local response for this, even for NonlinearPCM and SaLSA)
					const PCM& pcm = *((const PCM*)e.eVars.fluidSolver.get());
					ScalarField shapeSlab = getPlanarAvg(pcm.shape, iDir);
					for(int dir=0; dir<2; dir++)
						epsSlab[dir] += shapeSlab * (pcm.epsBulk - 1.); //fluid dielectric is isotropic
					kappaSqSlab += shapeSlab * pcm.k2factor;
					epsilonBulk = pcm.epsBulk;
					kappaSqBulk = pcm.k2factor;
				}
			}
			
			//Periodic potential and energy:
			ScalarFieldTilde dModel;
			bool periodicCoulomb = e.coulombParams.embedFluidMode || (!truncated); //if true, don't use slab truncation in periodic model
			Emodel = SlabPeriodicSolver(iDir, epsSlab[0], epsSlab[1], kappaSqSlab, periodicCoulomb).getEnergy(rhoModel, dModel);
			if(truncated) dModel = e.coulomb->embedShrink(dModel);
			Vmodel = I(dModel);
			
			//Isolated energy:
			std::shared_ptr<Coulomb> truncatedCoulomb;
			if(!truncated)
			{	//Still use truncation for the isolated case:
				CoulombParams truncatedParams = e.coulombParams;
				truncatedParams.geometry = geometry; //Slab
				truncatedParams.iDir = iDir;
				truncatedParams.embed = true;
				logSuspend();
				truncatedCoulomb = truncatedParams.createCoulomb(e.gInfo);
				logResume();
				#define EMBED_EXPAND(x, xBulk) \
				{	ScalarFieldTilde xTilde = J(x - xBulk); /*subtract bulk value*/ \
					xTilde = truncatedCoulomb->embedExpand(xTilde); /*switch to embedding grid*/ \
					x = xBulk + I(xTilde); \
				}
				for(int dir=0; dir<2; dir++) EMBED_EXPAND(epsSlab[dir], epsilonBulk)
				EMBED_EXPAND(kappaSqSlab, kappaSqBulk);
				#undef EMBED_EXPAND
			}
			CylindricalPoisson cp(iDir, epsSlab[0], epsSlab[1], kappaSqSlab);
			for(const Center& cdc: center)
			{	double zCenter = e.gInfo.R.column(iDir).length() * ws.restrict(cdc.pos - e.coulomb->xCenter)[iDir]; //Cartesian axial coordinate of center in embedding grid
				EmodelIsolated += cp.getEnergy(cdc.q, cdc.sigma, zCenter); //self energy of Gaussian (accounting for dielectric screening)
			}
			break;
		}
		default: die("\tCoulomb-interaction geometry must be either slab or periodic for charged-defect correction.\n");
	}
	logPrintf("\tEmodelIsolated: %.8lf\n", EmodelIsolated);
	logPrintf("\tEmodelPeriodic: %.8lf\n", Emodel);
	
	//Calculate alignment potential:
	ScalarField Varr[2] = { Vdft, Vmodel };
	vector3<> rCenter = e.gInfo.R*posMean;
	std::vector< std::vector<double> > hist = sphericalize(Varr, 2, 1., &rCenter);
	const std::vector<double>& rRadial = hist[0];
	const std::vector<double>& VdftRadial = hist[1];
	const std::vector<double>& VmodelRadial = hist[2];
	const std::vector<double>& wRadial = hist[3];
	//--- calculate alignment scalar:
	double deltaV = 0.; double wSum = 0.;
	for(size_t i=0; i<rRadial.size(); i++)
	{	double r = rRadial[i];
		double w = wRadial[i] * 0.5*erfc((rMin-r)/rSigma);
		deltaV += (VmodelRadial[i] - VdftRadial[i])*w;
		wSum += w;
	}
	deltaV /= wSum;
	logPrintf("\tDeltaV:         %.8lf (average over %lg grid points)\n", deltaV, wSum);
	logPrintf("\tNet correction: %.8lf (= EmodelIsolated - EmodelPeriodic + q deltaV)\n", EmodelIsolated - Emodel + qTot*deltaV);
	//--- write alignment data (spherical)
	string fname = e.dump.getFilename("chargedDefectDeltaV");
	logPrintf("\tWriting %s (spherically-averaged; plot to check DeltaV manually) ... ", fname.c_str()); logFlush();
	if(mpiUtil->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("\tError opening %s for writing.\n", fname.c_str())
		fprintf(fp, "#r DeltaV Vmodel Vdft weight\n");
		for(size_t i=0; i<rRadial.size(); i++)
			fprintf(fp, "%lf %le %le %le %le\n", rRadial[i], VmodelRadial[i]-VdftRadial[i],
				VmodelRadial[i], VdftRadial[i], (wRadial[i] * 0.5*erfc((rMin-rRadial[i])/rSigma)) / wSum);
		fclose(fp);
	}
	logPrintf("done.\n");
	//--- write alignment data (planar average)
	for(int iDir=0; iDir<3; iDir++)
	{	string fname = e.dump.getFilename(string("chargedDefectDeltaV") + "xyz"[iDir]);
		logPrintf("\tWriting %s (planarly-averaged normal to lattice direction# %d) ... ", fname.c_str(), iDir); logFlush();
		ScalarField Vavg[2]; double* VavgData[2];
		for(int k=0; k<2; k++)
		{	Vavg[k] = getPlanarAvg(Varr[k], iDir);
			VavgData[k] = Vavg[k]->data();
		}
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("\tError opening %s for writing.\n", fname.c_str())
			fprintf(fp, "#x[bohr] DeltaV Vmodel Vdft\n");
			vector3<int> iR;
			double h = e.gInfo.h[iDir].length();
			for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
			{	size_t i = e.gInfo.fullRindex(iR);
				fprintf(fp, "%lf %le %le %le\n", iR[iDir]*h,
					VavgData[1][i]-VavgData[0][i], VavgData[1][i], VavgData[0][i]);
			}
			fclose(fp);
		}
		logPrintf("done.\n");
	}
}
