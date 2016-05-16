#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/Dump_internal.h>
#include <core/WignerSeitz.h>
#include <core/Operators.h>
#include <core/ScalarFieldIO.h>
#include <core/LatticeUtils.h>
#include <core/Coulomb_internal.h>
#include <gsl/gsl_sf.h>

//-------------------------- Slab epsilon ----------------------------------

inline void planarAvg_sub(size_t iStart, size_t iStop, const vector3<int>& S, int iDir, complex* data)
{	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	THREAD_halfGspaceLoop( if(iG[jDir] || iG[kDir]) data[i] = 0.; )
}
void planarAvg(ScalarFieldTilde& X, int iDir)
{	threadLaunch(planarAvg_sub, X->gInfo.nG, X->gInfo.S, iDir, X->data());
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

void SlabEpsilon::dump(const Everything& e, ScalarField d_tot) const
{	string fname = e.dump.getFilename("slabEpsilon");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	//Read reference Dtot:
	ScalarField d_totRef(ScalarFieldData::alloc(e.gInfo));
	loadRawBinary(d_totRef, dtotFname.c_str());
	//Calculate inverse of epsilon:
	int iDir = e.coulombParams.iDir;
	vector3<> zHat = e.gInfo.R.column(iDir);
	zHat *= 1./zHat.length(); //unit vector along slab normal
	double dE = dot(zHat, e.coulombParams.Efield - Efield);
	if(!dE) die("\nThe applied electric fields in the reference and present calculations are equal.\n");
	//--- calculate field using central-difference derivative:
	ScalarFieldTilde tPlus(ScalarFieldTildeData::alloc(e.gInfo)), tMinus(ScalarFieldTildeData::alloc(e.gInfo));
	initTranslation(tPlus, e.gInfo.h[iDir]);
	initTranslation(tMinus, -e.gInfo.h[iDir]);
	double h = e.gInfo.h[iDir].length();
	ScalarFieldTilde epsInvTilde = (1./(dE * 2.*h)) * (tPlus - tMinus) * J(d_tot - d_totRef);
	planarAvg(epsInvTilde, iDir);
	//Fix values of epsilon near truncation boundary to 1:
	ScalarField epsInv = I(epsInvTilde);
	threadLaunch(fixBoundary_sub, e.gInfo.nr, e.gInfo.S, iDir, e.coulomb->ivCenter[iDir] + e.gInfo.S[iDir]/2, epsInv->data());
	//Apply smoothing:
	epsInv = I(gaussConvolve(J(epsInv), sigma));
	//Write file:
	FILE* fp = fopen(fname.c_str(), "w");
	fprintf(fp, "#distance[bohr]  epsilon\n");
	vector3<int> iR;
	double* epsInvData = epsInv->data();
	for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
		fprintf(fp, "%lf %lf\n", h*iR[iDir], 1./epsInvData[e.gInfo.fullRindex(iR)]);
	fclose(fp);
	logPrintf("done\n"); logFlush();
}

//-------------------------- Charged defects ----------------------------------

struct SlabPeriodicSolver : public LinearSolvable<ScalarFieldTilde>
{	int iDir; //truncated direction
	const ScalarField& epsilon;
	const GridInfo& gInfo;
	RealKernel Ksqrt, Kinv;
	const ScalarField epsInv;
	double K0; //G=0 component of slab kernel
	
	static inline void setKernels_sub(size_t iStart, size_t iStop, const GridInfo* gInfo, int iDir, double epsMean, double kRMS, double* Ksqrt, double* Kinv)
	{	const vector3<int>& S = gInfo->S;
		CoulombSlab_calc slabCalc(iDir, 0.5*gInfo->R.column(iDir).length());
		THREAD_halfGspaceLoop
		(	double K = slabCalc(iG, gInfo->GGT) / (4*M_PI);
			Kinv[i] = (fabs(K)>1e-12) ? 1./K : 0.;
			//Kinv[i] = gInfo->GGT.metric_length_squared(iG);
			Ksqrt[i] = (Kinv[i] || kRMS) ? 1./(epsMean*sqrt(fabs(Kinv[i]) + kRMS*kRMS)) : 0.;
		)
	}
	
	SlabPeriodicSolver(int iDir, const ScalarField& epsilon)
	: iDir(iDir), epsilon(epsilon), gInfo(epsilon->gInfo), Ksqrt(gInfo), Kinv(gInfo), epsInv(inv(epsilon))
	{
		threadLaunch(setKernels_sub, gInfo.nG, &gInfo, iDir, integral(epsilon)/gInfo.detR, 0., Ksqrt.data, Kinv.data);
		Ksqrt.set(); Kinv.set();
		nullToZero(state, gInfo);
	}
	
	ScalarFieldTilde hessian(const ScalarFieldTilde& phiTilde) const
	{	ScalarFieldTilde rhoTilde = -(Kinv * phiTilde); //vacuum term
		rhoTilde += divergence(J((epsilon-1.) * I(gradient(phiTilde))));  //dielectric term
		return (-1./(4*M_PI)) * rhoTilde;
	}
	
	ScalarFieldTilde precondition(const ScalarFieldTilde& rTilde) const
	{	return Ksqrt*(J(epsInv*I(Ksqrt*rTilde)));
	}

	double getEnergy(const ScalarFieldTilde& rho, ScalarFieldTilde& phi)
	{	MinimizeParams mp;
		mp.nDim = gInfo.nr;
		mp.nIterations = 20;
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
	matrix epsilonTilde, kappaSqTilde; diagMatrix G; //Fourier space operators
	
	CylindricalPoisson(int iDir, const ScalarField& epsilonSlab, const ScalarField& kappaSqSlab)
	{
		//Extract grid dimensions and bulk response:
		const GridInfo& gInfo = epsilonSlab->gInfo;
		NZ = gInfo.S[iDir];
		L = gInfo.R.column(iDir).length();
		vector3<int> iRbulk; iRbulk[iDir] = NZ/2;
		size_t iBulk = gInfo.fullRindex(iRbulk);
		epsilonBulk = epsilonSlab->data()[iBulk];
		kappaSqBulk = kappaSqSlab->data()[iBulk];
		
		//Initialize Fourier space operators along truncated direction:
		complexScalarFieldTilde epsilonSlabTilde = J(Complex(epsilonSlab - epsilonBulk));
		complexScalarFieldTilde kappaSqSlabTilde = J(Complex(kappaSqSlab - kappaSqBulk));
		std::vector<complex> epsilonDiagTilde(NZ), kappaSqDiagTilde(NZ);
		for(int iZ=0; iZ<NZ; iZ++)
		{	vector3<int> iG; iG[iDir]=iZ;
			size_t iSlab = gInfo.fullRindex(iG);
			epsilonDiagTilde[iZ] = epsilonSlabTilde->data()[iSlab];
			kappaSqDiagTilde[iZ] = kappaSqSlabTilde->data()[iSlab];
		}
		epsilonTilde.init(NZ,NZ);
		kappaSqTilde.init(NZ,NZ);
		for(int iZ=0; iZ<NZ; iZ++)
			for(int jZ=0; jZ<NZ; jZ++)
			{	int kZ = (jZ - iZ);
				if(kZ<0) kZ += NZ; //wrap kZ to [0,NZ)
				epsilonTilde.set(iZ,jZ, epsilonDiagTilde[kZ]);
				kappaSqTilde.set(iZ,jZ, kappaSqDiagTilde[kZ]);
			}
		G.resize(NZ);
		for(int iZ=0; iZ<NZ; iZ++)
			G[iZ] = (2*M_PI/L) * (iZ<NZ/2 ? iZ : iZ-NZ);
	}
    
	//Integrand
	double integrand(double k, double sigma, double z0) const
	{	//Initialize k-dependent operators and source terms:
		matrix OgTilde(NZ,1); //Gaussian source term
		diagMatrix KinvTilde(NZ); //Truncated Greens function
		double alpha = sqrt(k*k + kappaSqBulk/epsilonBulk);
		for(int iZ=0; iZ<NZ; iZ++)
		{	OgTilde.set(iZ,0, exp(-0.5*std::pow(G[iZ]*sigma,2))*cis(-G[iZ]*z0)); //note O cancels 1./L in derivation
			KinvTilde[iZ] = L*epsilonBulk*(alpha*alpha + G[iZ]*G[iZ])/(1. - exp(-0.5*alpha*L)*cos(0.5*L*G[iZ]));
		}
		matrix chiTilde = (1./L) * (G*epsilonTilde*G + kappaSqTilde + (k*k)*epsilonTilde);
		double Uk = trace(dagger(OgTilde) * inv(KinvTilde + chiTilde) * OgTilde).real();
		return k * exp(-std::pow(k*sigma,2)) * Uk;
	}
	struct IntegrandParams { double sigma, z0; const CylindricalPoisson* cp; };
	static double integrand_wrapper(double k, void* params) //wrapper for GSL integration routine
	{	const IntegrandParams& ip = *((const IntegrandParams*)params);
		return ip.cp->integrand(k, ip.sigma, ip.z0);
	}
	
	//Calculate self energy of Gaussian with norm q and width sigma centered at z0
	double getEnergy(double q, double sigma, double z0) const
	{	size_t wsSize = 1024;
		gsl_integration_workspace* ws = gsl_integration_workspace_alloc(wsSize);
		IntegrandParams ip = { sigma, z0, this };
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
	vector3<> pos0 = e.coulombParams.geometry==CoulombParams::Periodic ? center[0].pos : e.coulomb->xCenter;
	vector3<> posMean = pos0;
	for(const Center& cdc: center)
		posMean += (1./center.size()) * ws.restrict(cdc.pos - pos0);
	
	//Read reference Dtot and calculate electrostatic potential difference within DFT:
	ScalarField d_totRef(ScalarFieldData::alloc(e.gInfo));
	loadRawBinary(d_totRef, dtotFname.c_str());
	ScalarField Vdft = d_tot - d_totRef; //electrostatic potential of defect from DFT
	
	//Calculate isolated and periodic self-energy (and potential) of model charge
	ScalarField Vmodel; double Emodel=0., EmodelIsolated=0.;
	switch(e.coulombParams.geometry)
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
		{	int iDir = e.coulombParams.iDir;
			if(!e.coulombParams.embed)
				die("\tCoulomb truncation must be embedded for charged-defect correction in slab geometry.\n");
			rhoModel = e.coulomb->embedExpand(rhoModel); //switch to embedding grid
			//Create dielectric model for slab:
			ScalarField epsSlab; nullToZero(epsSlab, e.gInfo);
			double* epsSlabData = epsSlab->data();
			std::ifstream ifs(slabEpsFname.c_str());
			if(!ifs.is_open()) die("\tCould not open slab dielectric model file '%s' for reading.\n", slabEpsFname.c_str());
			string commentLine; getline(ifs, commentLine); //get and ignore comment line
			vector3<int> iR;
			double h = e.gInfo.h[iDir].length();
			for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
			{	double dExpected = h*iR[iDir];
				double d; ifs >> d >> epsSlabData[e.gInfo.fullRindex(iR)];
				if(fabs(d-dExpected) > symmThreshold)
					die("\tGeometry mismatch in '%s': expecting distance %lg on line %d; found %lg instead.\n",
						slabEpsFname.c_str(), dExpected, iR[iDir]+2, d);
			}
			ScalarFieldTilde epsSlabMinus1tilde = J(epsSlab) * (e.gInfo.nr / e.gInfo.S[iDir]); //multiply by number of points per plane (to account for average below)
			epsSlabMinus1tilde->setGzero(epsSlabMinus1tilde->getGzero() - 1.); //subtract 1
			planarAvg(epsSlabMinus1tilde, iDir); //now contains a planarly-uniform version of epsSlab-1
			epsSlab = 1. + I(e.coulomb->embedExpand(epsSlabMinus1tilde)); //switch to embedding grid (note embedding eps-1 (instead of eps) since it is zero in vacuum)
			
			//Periodic potential and energy:
			ScalarFieldTilde dModel;
			Emodel = SlabPeriodicSolver(iDir, epsSlab).getEnergy(rhoModel, dModel);
			Vmodel = I(e.coulomb->embedShrink(dModel));
			
			//Isolated energy:
			ScalarField kappaSqSlab; nullToZero(kappaSqSlab, epsSlab->gInfo);
			CylindricalPoisson cp(iDir, epsSlab, kappaSqSlab);
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
		{	ScalarFieldTilde Vtilde = J(Varr[k]);
			planarAvg(Vtilde, iDir);
			Vavg[k] = I(Vtilde);
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
					VavgData[1][i]-VavgData[0][i], VavgData[0][i], VavgData[1][i]);
			}
			fclose(fp);
		}
		logPrintf("done.\n");
	}
}
