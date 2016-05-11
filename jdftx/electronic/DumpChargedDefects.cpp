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
			Ksqrt[i] = (Kinv[i] || kRMS) ? 1./(epsMean*sqrt(Kinv[i] + kRMS*kRMS)) : 0.;
		)
	}
	
	SlabPeriodicSolver(int iDir, const ScalarField& epsilon)
	: iDir(iDir), epsilon(epsilon), gInfo(epsilon->gInfo), Ksqrt(gInfo), Kinv(gInfo), epsInv(inv(epsilon))
	{
		threadLaunch(setKernels_sub, gInfo.nG, &gInfo, iDir, integral(epsilon)/gInfo.detR, 0., Ksqrt.data, Kinv.data);
		K0 = (4*M_PI)/Kinv.data[0];
		Kinv.data[0] = 0.;
		Ksqrt.data[0] = 0.;
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
		phi->setGzero(K0 * rho->getGzero());
		return -0.5*dot(phi, O(hessian(phi))) + dot(phi, O(rho)); //first-order correct estimate of final answer
	}
};


struct CylindricalPoisson
{
	int NZ;
	diagMatrix z; //real-space grid
	matrix IZ, JZ, OZ; //axial operators (full-matrices)
	diagMatrix DZ, EZ; //axial operators (diagonal ones)

	CylindricalPoisson(int iDir, const ScalarField& epsSlab, double z0scale=1.)
	{
		//Determine grid mapping parameters:
		const GridInfo& gInfo = epsSlab->gInfo;
		double dzTyp = gInfo.h[iDir].length();
		double z0 = 0.25*gInfo.R.column(iDir).length() * z0scale;
		double dZtarget = 2.*dzTyp/(z0*(5.+sqrt(5.)));
		NZ = 2*ceil(1./dZtarget); //keep even (Z in (-1,1))
		logPrintf("\tSlabIsolated: Axial coordinate map with z0: %lg  NZ: %d\n", z0, NZ);
		
		//Initialize grid and spectral basis:
		diagMatrix Z(NZ); z.resize(NZ);
		for(int iZ=0; iZ<NZ; iZ++)
		{	Z[iZ] = -1. + (iZ+1)*(2./(NZ+1));
			z[iZ] = z0*Z[iZ]/(1.-Z[iZ]*Z[iZ]);
		}
		IZ.init(NZ, NZ);
		{	complex* IZdata = IZ.data();
			for(int alphaZ=0; alphaZ<NZ; alphaZ++)
				for(int iZ=0; iZ<NZ; iZ++)
					*(IZdata++) = sin((alphaZ+1)*(Z[iZ]+1.)*(0.5*M_PI));
		}
		
		//Calculate the analytical matrix element generators:
		diagMatrix OZtilde(2*NZ+1);
		for(int n=0; n<(2*NZ+1); n++)
			OZtilde[n] = -0.25* n*M_PI * gsl_sf_Si(n*M_PI);
		
		//Calculate the numerical matrix element generators:
		//--- initialize quadrature grid:
		const int Nquad = 11; //number of points in quadrature
		const int NZint = (NZ + 1) * Nquad;
		diagMatrix Zint(NZint), wZint(NZint);
		{	gsl_integration_glfixed_table* glTable = gsl_integration_glfixed_table_alloc(Nquad);
			double* ZintPtr = Zint.data();
			double* wZintPtr = wZint.data();
			double deltaZ = 2./(NZ+1);
			for(int iZ=0; iZ<=NZ; iZ++)
			{	double Zlo = -1. + deltaZ*iZ;
				double Zhi = -1. + deltaZ*(iZ+1);
				for(int iQuad=0; iQuad<Nquad; iQuad++)
					gsl_integration_glfixed_point(Zlo, Zhi, iQuad, ZintPtr++, wZintPtr++, glTable);
			}
			gsl_integration_glfixed_table_free(glTable);
		}
		//--- interpolate dielectric function to quadrature grid:
		diagMatrix epsInt(NZint); double epsBG;
		{	//Extract spline coefficients from epsSlab:
			int Sdir = gInfo.S[iDir];
			int Shlf = Sdir / 2;
			std::vector<double> epsSamples;
			vector3<int> iR;
			for(iR[iDir]=Sdir-Shlf; iR[iDir]<Sdir; iR[iDir]++)
				epsSamples.push_back(epsSlab->data()[gInfo.fullRindex(iR)]);
			for(iR[iDir]=0; iR[iDir]<=Shlf; iR[iDir]++)
				epsSamples.push_back(epsSlab->data()[gInfo.fullRindex(iR)]);
			std::vector<double> epsCoeff = QuinticSpline::getCoeff(epsSamples);
			double zStart = -Shlf * dzTyp;
			double zStop = +Shlf * dzTyp;
			double dzInv = 1./dzTyp;
			//Set and check background value:
			epsBG = epsSamples.front();
			assert(epsBG == epsSamples.back());
			//Evaluate on quadrature grid:
			FILE* fp = fopen("test_z.quad_eps", "w");
			for(int iZint=0; iZint<NZint; iZint++)
			{	double z = z0*Zint[iZint]/(1.-Zint[iZint]*Zint[iZint]);
				epsInt[iZint] = (z<=zStart || z>=zStop)
					? epsBG
					: QuinticSpline::value(epsCoeff.data(), (z-zStart)*dzInv);
				fprintf(fp, "%lg %lg %lg\n", Zint[iZint], wZint[iZint], epsInt[iZint]);
			}
			fclose(fp);
		}
		//--- compute integrands without the oscillatory factors:
		diagMatrix EZtildeInt(NZint), DZtildeInt(NZint);
		std::vector<diagMatrix> cosRefInt(2, diagMatrix(NZint));
		for(int iZint=0; iZint<NZint; iZint++)
		{	const double& Zcur = Zint[iZint];
			const double& epsCur = epsInt[iZint];
			double Zcomb = (1.+Zcur*Zcur) / std::pow(1.-Zcur*Zcur,2);
			EZtildeInt[iZint] = wZint[iZint] * 0.5 * Zcomb * (epsCur - epsBG);
			DZtildeInt[iZint] = wZint[iZint] * ((0.125*M_PI*M_PI) / Zcomb) * epsCur;
			cosRefInt[0][iZint] = 1.;
			cosRefInt[1][iZint] = cos(0.5*M_PI*(1.+Zcur));
		}
		//--- compute the integrals:
		diagMatrix EZtilde(2*NZ+1), DZtilde(2*NZ+1);
		diagMatrix cosInt(NZint);
		for(int n=0; n<(2*NZ+1); n++)
		{	//Set the oscillatory term:
			double kZ = (0.5*M_PI) * n;
			for(int iZint=0; iZint<NZint; iZint++)
				cosInt[iZint] = cos(kZ * (1.+Zint[iZint]));
			EZtilde[n] = dot(EZtildeInt, cosInt - cosRefInt[n%2]);
			DZtilde[n] = dot(DZtildeInt, cosInt);
		}
		
		//Calculate the matrix elements in spectral basis:
		matrix OZ0(NZ,NZ);
		matrix DZ0(NZ,NZ);
		matrix EZ0(NZ,NZ);
		for(int iZ=0; iZ<NZ; iZ++)
		for(int jZ=0; jZ<NZ; jZ++)
		{	int iDiff = abs(iZ-jZ);
			int iSum = iZ+jZ+2;
			double OZ0ij = iDiff%2==0 ? z0 * (OZtilde[iDiff] - OZtilde[iSum]) : 0.;
			OZ0.set(iZ,jZ, OZ0ij);
			EZ0.set(iZ,jZ, epsBG*OZ0ij + z0 * (EZtilde[iDiff] - EZtilde[iSum]));
			DZ0.set(iZ,jZ, ((iZ+1)*(jZ+1)/z0) * (DZtilde[iDiff] + DZtilde[iSum]));
		}
		
		//Switch to orthonormal basis that diagonalizes differential:
		//--- orthonormalize w.r.t EZ, and diagonalize DZ:
		matrix UZ = invsqrt(EZ0), DZevecs;
		(UZ * DZ0 * UZ).diagonalize(DZevecs, DZ); //DZ is now ready in diagonal form
		matrix VZ = UZ * DZevecs; //net transformation from spectral to new basis
		EZ = eye(NZ); //EZ is now identity by definition
		OZ = dagger(VZ) * OZ0 * VZ;
		IZ = IZ * VZ;
		JZ = inv(IZ);
	}
    
    //Return a normalized gaussian of width = sigma and integral = norm cenetred on the z-axis at zCenter:
    matrix gaussian(double norm, double sigma, double zCenter)
	{	double expFac = -0.5/(sigma*sigma);
		double normFac = norm/(sqrt(2.*M_PI)*sigma);
		matrix result(NZ,1);
		for(int iZ=0; iZ<NZ; iZ++)
			result.set(iZ,0, normFac*exp(expFac*std::pow(z[iZ]-zCenter,2)));
		return result;
	}
	
	//Calculate exp(x) * E1(x) handling overflow/underflow issues correctly
	double expE1prod(double x)
	{	if(x < 50.) return exp(x) * gsl_sf_expint_E1(x);
		else //Laurent series with ~ 1e-16 accuracy for x > 50.:
		{	double xInv = 1./x;
			double result = 1.;
			for(int n=20; n>0; n--)
				result = 1. - n * xInv * result;
			return xInv * result;
		}
	}

	
	double getEnergy(double q, double sigma, double z)
	{	matrix Orho = OZ * (JZ * gaussian(q, sigma, z));
		double U = 0.;
		for(int iZ=0; iZ<NZ; iZ++)
		{	double q = Orho(iZ,0).real(); //q_i from derivation, since we are in DZ-diagonal basis
			double betaSq = DZ[iZ]; //betaSq_i from derivation, since we are in DZ-diagonal basis
			U += 0.5*q*q * expE1prod(betaSq * sigma*sigma);
		}
		return U;
	}
};
//Refine CylindricalPoisson results with Richardson extrapolation:
struct SlabIsolated
{
	const vector3<> zScale; //scale factors
	vector3<std::shared_ptr<CylindricalPoisson> > cp;
	matrix3<> Uextrap;
	
	SlabIsolated(int iDir, const ScalarField& epsSlab) : zScale(0.8, 1.0, 1.25)
	{	//Initialize CylindircalPoisson instances:
		for(int i=0; i<3; i++)
			cp[i] = std::make_shared<CylindricalPoisson>(iDir, epsSlab, zScale[i]);
		//Initialize extrapolation matrix:
		matrix3<> scaleMat;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				scaleMat(i,j) = std::pow(zScale[i], -j);
		Uextrap = inv(scaleMat);
	}
	
	double getEnergy(double q, double sigma, double z)
	{	vector3<> U;
		for(int i=0; i<3; i++)
			U[i] = cp[i]->getEnergy(q, sigma, z);
		return (Uextrap * U)[0];
	}
};


//Get the averaged field in direction iDir between planes iCenter +/- iDist
double getEfield(ScalarField V, int iDir, int iCenter, int iDist)
{	const GridInfo& gInfo = V->gInfo;
	//Planarly average:
	ScalarFieldTilde Vtilde = J(V);
	planarAvg(Vtilde, iDir);
	ScalarField Vavg = I(Vtilde);
	//Extract field:
	assert(2*iDist < gInfo.S[iDir]);
	vector3<int> iRm, iRp;
	iRm[iDir] = positiveRemainder(iCenter - iDist, gInfo.S[iDir]);
	iRp[iDir] = positiveRemainder(iCenter + iDist, gInfo.S[iDir]);
	double Vm = Vavg->data()[gInfo.fullRindex(iRm)];
	double Vp = Vavg->data()[gInfo.fullRindex(iRp)];
	double dx = (2*iDist) * gInfo.h[iDir].length();
	return (Vm - Vp) / dx;
}

//Add electric field in direction iDir to given data
//x0 in lattice coordinates specifies where the added potential is zero
void addEfield_sub(size_t iStart, size_t iStop, const GridInfo* gInfo, int iDir, double Efield, vector3<> x0, double* V)
{	const vector3<int>& S = gInfo->S;
	double h = gInfo->h[iDir].length();
	double L = h * S[iDir];
	double r0 = x0[iDir] * L;
	THREAD_rLoop
	(	double r = h * iv[iDir] - r0;
		r -= L * floor(0.5 + r/L); //wrap to [-L/2,+L/2)
		V[i] += -Efield*r;
	)
}

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
			//--- fix up net electric field in output potential (increases accuracy of alignment):
			Vmodel = I(dModel);
			double EfieldDft = getEfield(Vdft, iDir, e.coulomb->ivCenter[iDir], e.gInfo.S[iDir]/2-2);
			double EfieldModel = getEfield(Vmodel, iDir, 0., e.gInfo.S[iDir]/2-2);
			vector3<> posMeanEmbed = Diag(e.coulomb->embedScale) * ws.restrict(posMean - e.coulomb->xCenter);
			threadLaunch(addEfield_sub, Vmodel->gInfo.nr, &(Vmodel->gInfo), iDir, EfieldDft-EfieldModel, posMeanEmbed, Vmodel->data());
			Vmodel = I(e.coulomb->embedShrink(J(Vmodel)));
			
			//Isolated energy:
			SlabIsolated si(iDir, epsSlab);
			for(const Center& cdc: center)
			{	double zCenter = e.gInfo.R.column(iDir).length() * ws.restrict(cdc.pos - e.coulomb->xCenter)[iDir]; //Cartesian axial coordinate of center in embedding grid
				EmodelIsolated += si.getEnergy(cdc.q, cdc.sigma, zCenter); //self energy of Gaussian (accounting for dielectric screening)
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
