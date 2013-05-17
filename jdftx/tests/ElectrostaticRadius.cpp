/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/Polarizability.h>
#include <electronic/VanDerWaals.h>
#include <electronic/operators.h>
#include <electronic/SpeciesInfo_internal.h>
#include <electronic/FluidSolver.h>
#include <commands/parser.h>
#include <core/DataMultiplet.h>
#include <core/DataIO.h>
#include <core/Util.cpp>
#include <core/WignerSeitz.h>
#include <gsl/gsl_sf_bessel.h>


struct MultipoleBasis
{
	const Basis& basis;
	vector3<> center, *centerPref; //center of cell in lattice coordinates
	double Rmax; //grid truncation radius
	double Gmax; //max G-vector in basis
	
	MultipoleBasis(const Basis& basis, vector3<> center, double Ecut)
	: basis(basis), center(center), Rmax(WignerSeitz(basis.gInfo->R).inRadius()), Gmax(sqrt(2*Ecut))
	{	
		#ifdef GPU_ENABLED
		cudaMalloc(&centerPref, sizeof(vector3<>));
		cudaMemcpy(centerPref, &center, sizeof(vector3<>), cudaMemcpyHostToDevice);
		gpuErrorCheck();
		#else
		centerPref = &center;
		#endif
	}
	
	~MultipoleBasis()
	{
		#ifdef GPU_ENABLED
		cudaFree(centerPref);
		#endif
	}
	
	//Create the projectors for each G (> 0) in Garr at given l,m:
	ColumnBundle getProjector(int l, int m, const std::vector<double>& Garr)
	{	ColumnBundle proj(Garr.size(), basis.nbasis, &basis, 0, isGpuEnabled());
		complex* projData = proj.dataPref();
		//Initialize radial grids:
		const double dG(0.01 * (2*M_PI)/Rmax);
		const int nG = int(ceil(Gmax/dG))+5;
		const double sigma = 0.5, eta = sqrt(0.5)/sigma; //smooth dropoff at end of grid
		const double rStart = 1e-3, rStop = Rmax + 10*sigma, dlogr = 0.01, rRatio = exp(dlogr);
		std::vector<double> rArr; for(double r=rStart; r<=rStop; r*=rRatio) rArr.push_back(r);
		RadialFunctionR jl(rArr, dlogr);
		for(unsigned i=0; i<Garr.size(); i++)
		{	const double& G = Garr[i];
			assert(G > 0);
			//Create radial part in real space:
			for(unsigned j=0; j<rArr.size(); j++)
			{	const double& r = rArr[j];
				jl.f[j] = bessel_jl(l,G*r) * 0.5*erfc(eta*(r-Rmax));
			}
			//Bessel transform to reciprocal space:
			RadialFunctionG jlTilde;
			jl.transform(l, dG, nG, jlTilde);
			//Store onto projector:
			callPref(Vnl)(basis.nbasis, 0, 1, l, m, vector3<>(0,0,0), basis.iGarrPref, basis.gInfo->G,
				centerPref, jlTilde, projData+proj.index(i,0), false, vector3<complex*>());
			jlTilde.free();
		}
		proj *= (sqrt(4*M_PI) * cis(-l*M_PI/2)); //makes stuff real for real wavefunctions and set normalization
		return proj;
	}
};

inline void rArrSet(size_t iStart, size_t iStop, vector3<int> S, matrix3<> R, vector3<> center, vector3<double*> rArr)
{	matrix3<> invS = inv(matrix3<>(Diag(S)));
	THREAD_rLoop
	(	vector3<> x = invS * iv;
		for(int j=0; j<3; j++)
			x[j] -= floor(x[j]-center[j]+0.5);
		storeVector(R * (x - center), rArr, i);
	)
}

void reduceGplanar(const DataGptr& x, matrix& xRed)
{	const GridInfo& gInfo = x->gInfo;
	int nG = xRed.nRows()/2;
	assert(xRed.nCols() == 1);
	assert(gInfo.S[0] == 1);
	assert(gInfo.S[1] == 1);
	assert(gInfo.S[2] > nG);
	callPref(eblas_copy)(xRed.dataPref(), x->dataPref()+1, nG);
	callPref(eblas_copy)(xRed.dataPref()+nG, x->dataPref()+1, nG);
	callPref(eblas_dscal)(nG, -1., ((double*)(xRed.dataPref()+nG))+1, 2); //negate the imaginary parts (complex conjugate if inversion symmetry employed)
}

DataRptr IexpandGplanar(const matrix& xRed, const GridInfo& gInfo, double sigma=0.)
{	int nG = xRed.nRows()/2;
	assert(xRed.nCols() == 1);
	assert(gInfo.S[0] == 1);
	assert(gInfo.S[1] == 1);
	assert(gInfo.S[2] > nG);
	DataGptr x; nullToZero(x, gInfo);
	callPref(eblas_copy)(x->dataPref()+1, xRed.dataPref(), nG);
	if(sigma) x = gaussConvolve(x, sigma);
	return I(x);
}

int main(int argc, char** argv)
{	if(argc < 2) { logPrintf("Usage: SphericalChi <jdftx-input-file> [<r0>=0] [<r1>=0] [<r2>=0]\n"); return 1; }
	initSystem(argc, argv);
	
	//Read the input file:
	Everything e;
	parse(argv[1], e);
	//Skip reading the wavefunctions and initialzing empty states to save time:
	assert(e.dump.polarizability);
	Polarizability& pol = *e.dump.polarizability;
	int nV = e.eInfo.nElectrons/2;
	int nC = e.eInfo.nBands - nV;
	if(pol.Ecut<=0.) pol.Ecut = 4.*e.cntrl.Ecut;
	if(pol.nEigs<=0) pol.nEigs = nV * nC;
	e.eInfo.nBands = nV;
	e.eVars.wfnsFilename.clear();
	//Perform JDFTx initialization:
	e.setup();
	
	//read molecule center fromcommandline:
	vector3<> center; //molecule center in cartesian coordinates
	for(int j=0; j<3; j++) if(argc>j+2) center[j] = atof(argv[j+2]);
	logPrintf("\nMolecule center at "); center.print(globalLog, " %lg ");
	center = e.gInfo.invR * center; //convert to lattice coordinates
	
	//Setup the basis and read in polarizability:
	logPrintf("Setting up reduced basis at Ecut=%lg: ", pol.Ecut);
	Basis basis; basis.setup(e.gInfo, e.iInfo, pol.Ecut, pol.dk);
	logPrintf("Reading polarizability eigensystem\n");
	ColumnBundle V(pol.nEigs, basis.nbasis, &basis);
	matrix Xext(pol.nEigs, pol.nEigs);
	((ManagedMemory&)V).read(e.dump.getFilename("pol_basis").c_str());
	Xext.read(e.dump.getFilename("pol_Xext").c_str());
	
	//Diagonalize Xext:
	logPrintf("Diagonalizing polarizability\n");
	{	matrix Xext_evecs; diagMatrix Xext_eigs;
		Xext.diagonalize(Xext_evecs, Xext_eigs);
		diagMatrix Xext_eigsSqrt = Xext_eigs; for(double& x: Xext_eigsSqrt) { assert(x<0); x = sqrt(-x); }
		V = V * (Xext_evecs * Xext_eigsSqrt);
	}
	
	//Compute net charge density (with smoothed nuclei ala nonlocal PCM):
	logPrintf("Adding rotational modes\n");
	DataGptr rho = J(e.eVars.get_nTot());
	logSuspend(); VanDerWaals vdW(e); logResume();
	for(const auto& sp: e.iInfo.species)
	{	DataGptr SG; nullToZero(SG, e.gInfo);
		callPref(getSG)(e.gInfo.S, sp->atpos.size(), sp->atposPref, 1./e.gInfo.detR, SG->dataPref()); //get structure factor for current species
		rho -= sp->Z * gaussConvolve(SG, vdW.getParams(sp->atomicNumber).R0/6.);
	}
	//Replace the last (least significant) polarizability eigencomponent with the rotational one:
	V.setColumn(V.nCols()-1, sqrt(1./e.eVars.fluidParams.T)*Complex(rho));
	
	//Get liquid properties:
	if(e.eVars.fluidParams.solvents.size()!=1) die("Fluid must include exactly one solvent component.\n");
	double Nbulk = e.eVars.fluidParams.solvents[0]->Nbulk;
	double epsBulk = e.eVars.fluidSolver->epsBulk;
	DataRptrVec rArr(e.gInfo); threadLaunch(rArrSet, e.gInfo.nr, e.gInfo.S, e.gInfo.R, center, rArr.data()); //Create array of r in real space
	ColumnBundle OJr(3, basis.nbasis, &basis); for(int k=0; k<3; k++) OJr.setColumn(k, O(J(Complex(rArr[k])))); //Convert to a columnbundle projector
	//Compute mean field dielectric constant
	double epsMF = 1. + 4*M_PI * Nbulk * pow(nrm2(OJr ^ V),2)/3.;
	double Aeps = 1. + 1./(epsBulk-1.) - 1./(epsMF-1.);
	logPrintf("Nbulk: %lg   epsBulk: %lg   epsMF: %lg   Aeps: %lg\n", Nbulk, epsBulk, epsMF, Aeps);

	//Create capacitor geometry (1D planar):
	double L = 384;
	double dG = (2*M_PI)/L;
	//Reduced 'sphere':
	double Gmax = sqrt(2*pol.Ecut);
	int nG = int(ceil(Gmax/dG));
	std::vector<double> Garr(nG);
	for(int iG=0; iG<nG; iG++) Garr[iG] = (iG+1)*dG;
	//Density grid:
	GridInfo gPlanar;
	gPlanar.R = Diag(vector3<>(1., 1., L));
	gPlanar.S = vector3<int>(1, 1, 3840);
	gPlanar.initialize();
	double dz = L / gPlanar.nr;
	const double xHlf = 0.25 - (15.+45./Gmax)/L; //half-length of liquid in lattice coordinates
	logPrintf("Liquid boundary at z = %lg\n", xHlf * L);
	
	//Compute the susceptibility matrix (in 1D GG'):
	matrix chi = zeroes(2*nG,2*nG);
	logSuspend(); MultipoleBasis mb(basis, center, pol.Ecut); logResume();
	for(int l=0; l<=3; l++)
		for(int m=-l; m<=+l; m++)
		{	logPrintf("Adding contributions to chi at l=%d, m=%+d.\n", l, m); logFlush();
			matrix Vlm = mb.getProjector(l,m,Garr) ^ V;
			matrix chiHlf = Vlm * dagger(Vlm); //a GG' matrix in the G>0 G'>0 quadrant
			//Add in the other quadrants:
			int lSign = (l%2==0 ? +1 : -1);
			matrix chiSub(2*nG,2*nG);
			chiSub.set(0,nG, 0,nG, chiHlf);
			chiSub.set(nG,2*nG, 0,nG, lSign*chiHlf);
			chiSub.set(0,nG, nG,2*nG, lSign*dagger(chiHlf));
			chiSub.set(nG,2*nG, nG,2*nG, chiHlf);
			chi += chiSub;
		}
	
	//Multiply by the shape function:
	logPrintf("Applying shape function.\n");
	complex* chiData = chi.data();
	for(int i=0; i<2*nG; i++)
	{	int iG = (i<nG ? i+1 : -(i-nG+1));
		for(int j=0; j<2*nG; j++)
		{	int jG = (j<nG ? j+1 : -(j-nG+1));
			double deltaGL = (iG - jG) * (2.*M_PI);
			double sTilde = 2.* (deltaGL ? sin(xHlf*deltaGL)/deltaGL : xHlf);
			chiData[chi.index(i,j)] *= (-sTilde * Nbulk);
		}
	}
	
	//Create truncated Coulomb kernel:
	diagMatrix K(2*nG);
	for(int i=0; i<2*nG; i++)
	{	double G = dG * (i<nG ? i+1 : -(i-nG+1));
		K[i] = 4*M_PI*(1.-cos(0.5*G*L))/(G*G);
	}
	
	//Create test charge density:
	DataRptr rhoExt; nullToZero(rhoExt, gPlanar);
	rhoExt->data()[gPlanar.S[2]/4] = 1./gPlanar.dV;
	rhoExt->data()[(3*gPlanar.S[2])/4] = -1./gPlanar.dV;
	//DataGptr rhoExtTilde = gaussConvolve(J(rhoExt), 9./Gmax); //make resolvable on reduced grid
	//rhoExt = I(rhoExtTilde);
	
	matrix rhoExtMat(2*nG, 1);
	reduceGplanar(J(rhoExt), rhoExtMat);
	matrix phiExtMat = K * rhoExtMat;

	matrix rhoBoundMat = dagger_symmetrize(chi * inv(eye(2*nG) - Aeps*K*chi)) * phiExtMat;
	matrix phiTotMat = K * rhoBoundMat + phiExtMat;
	
	double sigma = 3./Gmax;
	DataRptr phiExt = IexpandGplanar(phiExtMat, gPlanar, sigma);
	DataRptr rhoBound = IexpandGplanar(rhoBoundMat, gPlanar, sigma);
	DataRptr phiTot = IexpandGplanar(phiTotMat, gPlanar, sigma);
	
	FILE* fp = fopen(e.dump.getFilename("elRadiusPlot").c_str(), "w");
	double *rhoExtData = rhoExt->data(), *phiExtData = phiExt->data();
	double *rhoBoundData = rhoBound->data(), *phiTotData = phiTot->data();
	fprintf(fp, "#z[bohr]\trhoExt\tphiExt\trhoBound\tphiTot\n");
	for(int i=0; i<gPlanar.nr; i++)
		fprintf(fp, "%lg\t%le\t%le\t%le\t%le\n", i*dz, rhoExtData[i], phiExtData[i], rhoBoundData[i], phiTotData[i]);
	fclose(fp);
	
	double rhoSum = 0., rhoSumZ = 0.;
	for(int i=0; i<gPlanar.nr/4; i++)
	{	rhoSum += rhoBoundData[i];
		rhoSumZ += rhoBoundData[i] * (i*dz - xHlf*L);
	}
	logPrintf("Electrostatic radius = %lg bohrs.\n", rhoSumZ / rhoSum);
	
	finalizeSystem();
	return 0;
}
