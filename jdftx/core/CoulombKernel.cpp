/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <core/CoulombKernel.h>
#include <core/Coulomb_internal.h>
#include <core/LatticeUtils.h>
#include <core/LoopMacros.h>
#include <core/ManagedMemory.h>
#include <core/Thread.h>
#include <cfloat>

const double CoulombKernel::nSigmasPerWidth = 1.+sqrt(-2.*log(DBL_EPSILON)); //gaussian negligible at double precision (+1 sigma for safety)

CoulombKernel::CoulombKernel(const matrix3<> R, const vector3<int> S, const vector3<bool> isTruncated, double omega)
: R(R), S(S), isTruncated(isTruncated), omega(omega)
{
}


void CoulombKernel::compute(double* data, const WignerSeitz& ws, symmetricMatrix3<>* data_RRT) const
{	//Count number of truncated directions:
	int nTruncated = 0;
	for(int k=0; k<3; k++) if(isTruncated[k]) nTruncated++;
	//Call appropriate routine:
	switch(nTruncated)
	{	case 2: computeWire(data, ws, data_RRT); break;
		case 3: computeIsolated(data, ws, data_RRT); break;
		default: assert(!"Invalid truncated direction count");
	}
}

//! Compute erfc(omega r)/r - erfc(a r)/r
inline double screenedPotential(double r, double a, double omega)
{	double result = a * erf_by_x(a * r);
	if(omega) result -= omega * erf_by_x(omega * r);
	return result;
}

//! Compute (1/r) d/dr of screenedPotential()
inline double screenedPotentialPrime_by_r(double r, double a, double omega)
{	double result = std::pow(a,3) * erf_by_xPrime_by_x(a * r);
	if(omega) result -= std::pow(omega,3) * erf_by_xPrime_by_x(omega * r);
	return result;
}

//--------- Fully truncated Coulomb Kernel ---------

namespace CoulombKernelIsolated
{
	//Initialize the long range part of the kernel in real space with the minimum image convention:
	inline void realSpace_thread(size_t iStart, size_t iStop, vector3<int> Sdense, matrix3<> R,
		double* data, const WignerSeitz* ws, double sigma, double omega,
		bool computeStress, size_t arrayStride)
	{
		vector3<> invSdense; for(int k=0; k<3; k++) invSdense[k] = 1./Sdense[k];
		double dV = fabs(det(R)) * (invSdense[0]*invSdense[1]*invSdense[2]); //integration factor
		matrix3<> RTR = (~R)*R; //metric
		double a = sqrt(0.5)/sigma;
		int nArr = computeStress ? 7 : 1;
		vector3<int> S = Sdense; S[2] = 2*(Sdense[2]/2+1); //padding for in-place r2c transform
		THREAD_rLoop
		(	if(iv[2]<Sdense[2]) 
			{	vector3<> x; for(int k=0; k<3; k++) x[k] = invSdense[k] * iv[k]; //lattice coordinates
				x = ws->reduce(x); //minimum image convention
				double r = sqrt(RTR.metric_length_squared(x)); //minimum image distance
				data[i] = dV * screenedPotential(r, a, omega);
				if(computeStress)
				{	vector3<> rVec = R * x;
					symmetricMatrix3<> V_RRT = (dV * screenedPotentialPrime_by_r(r, a, omega)) * outer(rVec);
					const double* V_RRTcomponent = (const double*)&V_RRT;
					for(int iComp=0; iComp<6; iComp++)
						data[i+(iComp+1)*arrayStride] = V_RRTcomponent[iComp]
							+ (iComp<3 ? data[i] : 0.); //contribution due to det(R) in dV
				}
			}
			else for(int iArr=0; iArr<nArr; iArr++) data[i+iArr*arrayStride] = 0.; //padded points
		)
	}

	//Add short-ranged parts in reciprocal space (and optionally downsample)
	inline void recipSpace_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
		const complex* in, const vector3<int>& Sdense, double* out, double sigma,
		const matrix3<>& G, symmetricMatrix3<>* out_RRT, size_t arrayStride)
	{
		vector3<int> pitchDense;
		pitchDense[2] = 1;
		pitchDense[1] = pitchDense[2] * (1+Sdense[2]/2);
		pitchDense[0] = pitchDense[1] * Sdense[1];
		double hlfSigmaSq = 0.5*sigma*sigma;
		THREAD_halfGspaceLoop
		(	//Find corresponding point in dense cell fourier transform:
			int iDense = 0;
			for(int k=0; k<3; k++)
				iDense += pitchDense[k] * (iG[k]<0 ? iG[k]+Sdense[k] : iG[k]);
			//Store in smaller cell, along with short-ranged part:
			double Gsq = GGT.metric_length_squared(iG);
			out[i] = in[iDense].real() + (4*M_PI) * (Gsq ? (1.-exp(-hlfSigmaSq*Gsq))/Gsq : hlfSigmaSq);
			//Optional lattice derivative:
			if(out_RRT)
			{	symmetricMatrix3<>& V_RRT = out_RRT[i];
				//Set the reciprocal-space contribution:
				vector3<> Gcart = iG * G;
				double minus_Vprime_by_G = Gsq ? (8*M_PI) * (1.-exp(-hlfSigmaSq*Gsq)*(1.+hlfSigmaSq*Gsq))/(Gsq*Gsq) : 0.;
				V_RRT = (minus_Vprime_by_G) * outer(Gcart);
				//Add the part computed from real space:
				double* V_RRTcomponent = (double*)&V_RRT;
				for(int iComp=0; iComp<6; iComp++)
					V_RRTcomponent[iComp] += in[iDense+(iComp+1)*arrayStride].real();
			}
		)
	}
}

void CoulombKernel::computeIsolated(double* data, const WignerSeitz& ws, symmetricMatrix3<>* data_RRT) const
{
	for(int k=0; k<3; k++) assert(isTruncated[k]); //Make sure all directions are truncated
	double sigma = ws.inRadius() / nSigmasPerWidth;
	logPrintf("Gaussian width for range separation: %lg bohrs.\n", sigma);
	
	//Set up dense integration grids:
	logPrintf("FFT grid for long-range part: ");
	double Gnyq = nSigmasPerWidth * std::max(1./sigma, omega*sqrt(2.)); //lower bound on Nyquist frequency
	vector3<int> Sdense; //dense fft sample count
	for(int k=0; k<3; k++)
	{	Sdense[k] = std::max(S[k], 2*int(ceil(Gnyq * R.column(k).length() / (2*M_PI))));
		while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
	}
	logPrintf("[ %d %d %d ].\n", Sdense[0], Sdense[1], Sdense[2]);
	size_t nG = S[0] * (S[1] * size_t(1 + S[2]/2)); //number of symmetry reduced reciprocal lattice vectors
	size_t nGdense = Sdense[0] * (Sdense[1] * size_t(1+Sdense[2]/2));
	size_t nArrays = data_RRT ? 7 : 1;
	size_t arrayStride = nGdense; //for the complex arrays
	size_t arrayStrideReal = 2*nGdense; //for the real arrays
	
	//Plan Fourier transforms:
	ManagedArray<complex> dense; dense.init(nGdense*nArrays);
	complex* denseArr = dense.data();
	double* denseRealArr = (double*)denseArr;
	logPrintf("Planning fourier transform ... "); logFlush();
	fftw_plan_with_nthreads(nProcsAvailable);
	fftw_plan fftPlanR2C = fftw_plan_dft_r2c_3d(Sdense[0], Sdense[1], Sdense[2], denseRealArr, (fftw_complex*)denseArr, FFTW_ESTIMATE);
	logPrintf("Done.\n");
	
	//Long-range part in real space
	logPrintf("Computing truncated long-range part in real space ... "); logFlush();
	threadLaunch(CoulombKernelIsolated::realSpace_thread, 2*nGdense, Sdense, R, denseRealArr, &ws, sigma, omega, data_RRT, arrayStrideReal);
	logPrintf("Done.\n");
	
	//Add short-ranged part in reciprocal space (and down-sample if required):
	logPrintf("Adding short-range part in reciprocal space ... "); logFlush();
	fftw_execute(fftPlanR2C);
	if(nArrays > 1)
	{	for(size_t iArray=1; iArray<nArrays; iArray++)
		{	double* arrayPtr = denseRealArr + iArray*arrayStrideReal;
			fftw_execute_dft_r2c(fftPlanR2C, arrayPtr, (fftw_complex*)arrayPtr);
		}
	}
	matrix3<> G = (2.*M_PI) * inv(R);
	matrix3<> GGT = G * (~G);
	threadLaunch(CoulombKernelIsolated::recipSpace_thread, nG, S, GGT, denseArr, Sdense, data, sigma, G, data_RRT, arrayStride);
	fftw_destroy_plan(fftPlanR2C);
	logPrintf("Done.\n");
}

//--------- Wire geometry (two directions truncated) coulomb kernel ---------

//Threaded computation for wire-geometry kernel (threaded over planes)
struct CoulombKernelWire
{	
	//Data arrays:
	ManagedArray<complex> dense; //2D fft array on dense grid
	double* Vc; //output kernel (3D fftw c2r layout)
	symmetricMatrix3<>* Vc_RRT; //lattice derivative of Vc
	fftw_plan fftPlanR2C; //2D fftw r2c plan
	
	//Geometry:
	matrix3<> R; vector3<int> S, Sdense; //lattice vector and sample counts
	int iDir, jDir, kDir; //iDir is the untruncated direction
	const WignerSeitz *ws; //Wigner-Seitz cell
	double sigma, omega;
	size_t arrayStride, arrayStrideReal;
	
	void computePlane(int iPlane)
	{
		//Geometry
		double L = R.column(iDir).length(); //Periodicity in untruncated direction
		double kCur = iPlane * (2.*M_PI)/L; //Wave vector of current plane along untruncated direction
		double rhoMax = ws->circumRadius(iDir);
		symmetricMatrix3<> axisHatOuter = outer(R.column(iDir)*(1./L)); //outer product of unit vector along periodic direction
		
		//Long range part in real space (on optionally denser grid):
		int jPitchDense = 2*(1+Sdense[kDir]/2);
		double invSjDense = 1./Sdense[jDir], invSkDense = 1./Sdense[kDir];
		double dA = fabs(det(R)) * (invSjDense*invSkDense) / L; //2D integration factor
		matrix3<> RTR = (~R)*R;
		Cbar_k_sigma cbar_k_sigma(kCur, sigma, rhoMax), *cbar_k_screen=0; //Look-up table for convolved cylindrical potential
		if(omega) cbar_k_screen = new Cbar_k_sigma(kCur, sqrt(0.5)/omega, rhoMax); //Look-up table for screened cylindrical potential
		Cbar_k_sigma *minus_cbar_k_sigma_k=0, *minus_cbar_k_screen_k=0; //Look-up tables for lattice derivatives
		if(Vc_RRT && iPlane)
		{	minus_cbar_k_sigma_k = new Cbar_k_sigma(kCur, sigma, rhoMax, 1., true);
			if(omega) minus_cbar_k_screen_k = new Cbar_k_sigma(kCur, sqrt(0.5)/omega, rhoMax, 1., true);
		}
		complex* denseArr = dense.data();
		double* denseRealArr = (double*)denseArr; //in-place transform
		for(int ij=0; ij<Sdense[jDir]; ij++)
			for(int ik=0; ik<Sdense[kDir]; ik++)
			{	vector3<> x; //point in lattice coordinates
				x[iDir] = 0.;
				x[jDir] = invSjDense * ij;
				x[kDir] = invSkDense * ik;
				x = ws->reduce(x);
				double rho = sqrt(RTR.metric_length_squared(x)); //minimum image distance in 2D
				int iDense = ik + jPitchDense * ij; //index into dense array (in the fftw in-place r2c layout)
				double term = dA * (cbar_k_sigma.value(rho) - (omega ? cbar_k_screen->value(rho) : 0.));
				denseRealArr[iDense] = term;
				if(Vc_RRT)
				{	double term_rho_by_rho=0., minus_term_k=0.;
					if(rho) term_rho_by_rho = dA * (cbar_k_sigma.deriv(rho) - (omega ? cbar_k_screen->deriv(rho) : 0.)) / rho;
					if(iPlane) minus_term_k = dA * (minus_cbar_k_sigma_k->value(rho) - (omega ? minus_cbar_k_screen_k->value(rho) : 0.));
					double curV_RRTzz = minus_term_k * kCur - term; //second term subtracts L contribution in dV propagation added below
					symmetricMatrix3<> curV_RRT = term_rho_by_rho * outer(R * x);
					const double* curV_RRTcomponent = (const double*)&curV_RRT;
					const double* axisHatOuterComponent = (const double*)&axisHatOuter;
					for(int iComp=0; iComp<6; iComp++)
						denseRealArr[iDense+(iComp+1)*arrayStrideReal]
							= curV_RRTcomponent[iComp]
							+ curV_RRTzz * axisHatOuterComponent[iComp]
							+ (iComp<3 ? term : 0.); //contribution via dV (writing dA as dV/L for simplicity)
				}
			}
		if(omega) delete cbar_k_screen;
		if(Vc_RRT && iPlane)
		{	delete minus_cbar_k_sigma_k;
			if(omega) delete minus_cbar_k_screen_k;
		}
		//Fourier transforms:
		fftw_execute_dft_r2c(fftPlanR2C, denseRealArr, (fftw_complex*)denseArr);
		if(Vc_RRT)
		{	for(size_t iArray=1; iArray<7; iArray++)
			{	double* arrayPtr = denseRealArr + iArray*arrayStrideReal;
				fftw_execute_dft_r2c(fftPlanR2C, arrayPtr, (fftw_complex*)arrayPtr);
			}
		}
		
		//Add analytic short-ranged parts in fourier space (and down-sample to final resolution):
		jPitchDense = 1+Sdense[kDir]/2;
		vector3<int> pitch;
		pitch[2] = 1;
		pitch[1] = pitch[2] * (1 + S[2]/2);
		pitch[0] = pitch[1] * S[1];
		matrix3<> G = (2.*M_PI) * inv(R);
		matrix3<> GGT = G * (~G);
		double hlfSigmaSq = 0.5*sigma*sigma;
		vector3<int> iG; iG[iDir] = iPlane;
		for(iG[jDir]=1-S[jDir]/2; iG[jDir]<=S[jDir]/2; iG[jDir]++)
			for(iG[kDir]=0; iG[kDir]<=S[kDir]/2; iG[kDir]++)
			{	//Collect the data from the dense grid transform:
				int iDense = iG[kDir] + jPitchDense*(iG[jDir]<0 ? iG[jDir]+Sdense[jDir] : iG[jDir]);
				double curV = denseArr[iDense].real();
				//Add the analytical short-ranged part:
				double Gsq = GGT.metric_length_squared(iG);
				curV += (4*M_PI) * (Gsq ? (1.-exp(-hlfSigmaSq*Gsq))/Gsq : hlfSigmaSq);
				//Compute corresponding lattice derivative if needed:
				symmetricMatrix3<> curV_RRT;
				if(Vc_RRT)
				{	//Set the reciprocal-space contribution:
					double minus_Vprime_by_G = Gsq ? (8*M_PI) * (1.-exp(-hlfSigmaSq*Gsq)*(1.+hlfSigmaSq*Gsq))/(Gsq*Gsq) : 0.;
					curV_RRT = (minus_Vprime_by_G) * outer(iG * G);
					//Add the part computed from real space:
					double* curV_RRTcomponent = (double*)&curV_RRT;
					for(int iComp=0; iComp<6; iComp++)
						curV_RRTcomponent[iComp] += denseArr[iDense+(iComp+1)*arrayStride].real();
				}
				//Save to the appropriate locations in Vc (4 sign combinations):
				for(int si=0; si<2; si++)
				{	for(int sk=0; sk<2; sk++)
					{	vector3<int> iv = iG;
						if(si && iv[iDir]) { iv[iDir] = S[iDir] - iv[iDir]; }
						if(sk && iv[kDir]) { iv[kDir] = S[kDir] - iv[kDir]; iv[jDir] = -iv[jDir]; }
						if(iv[jDir] < 0) iv[jDir] += S[jDir];
						if(iv[2] <= S[2]/2)
						{	size_t iOut = dot(iv, pitch);
							Vc[iOut] = curV;
							if(Vc_RRT) Vc_RRT[iOut] = curV_RRT;
						}
					}
				}
			}
	}
	
	static void thread(int iThread, int nThreads, CoulombKernelWire* ckwArr,
		int nPlanes, int* nPlanesDone, std::mutex* m)
	{
		while(true)
		{	//Get next available job:
			m->lock();
			int iPlane = (*nPlanesDone)++;
			m->unlock();
			if(iPlane >= nPlanes)
				break; //job queue empty
			//Perform job:
			ckwArr[iThread].computePlane(iPlane);
		}
	}
};

void CoulombKernel::computeWire(double* data, const WignerSeitz& ws, symmetricMatrix3<>* data_RRT) const
{	
	//Find axis and check geometry:
	int iDir = -1;
	for(int k=0; k<3; k++)
		if(!isTruncated[k])
		{	assert(iDir < 0); //there should be at most one untruncated direction
			iDir = k; //axis is along untruncated direction
		}
	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	assert(WignerSeitz::isOrthogonal(R.column(iDir), R.column(jDir)));
	assert(WignerSeitz::isOrthogonal(R.column(iDir), R.column(kDir)));
	double sigma = ws.inRadius(iDir) / nSigmasPerWidth;
	logPrintf("Gaussian width for range separation: %lg bohrs.\n", sigma);
	
	//Set up dense integration grids:
	logPrintf("FFT grid for long-range part: ");
	double Gnyq = nSigmasPerWidth * std::max(1./sigma, omega*sqrt(2.)); //lower bound on Nyquist frequency
	vector3<int> Sdense; //dense fft sample count
	for(int k=0; k<3; k++)
		if(k == iDir) Sdense[k] = S[k];
		else
		{	Sdense[k] = std::max(S[k], 2*int(ceil(Gnyq * R.column(k).length() / (2*M_PI))));
			while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
		}
	string dirName(3,'0'); dirName[iDir]='1';
	logPrintf("%d x %d perpendicular to %s direction.\n", Sdense[jDir], Sdense[kDir], dirName.c_str());
	size_t nGdensePlanar = Sdense[jDir] * size_t(1+Sdense[kDir]/2); //Number of wave vectors per densely sampled plane
	size_t nArrays = data_RRT ? 7 : 1;

	//Plan Fourier transforms:
	logPrintf("Planning fourier transform ... "); logFlush();
	fftw_plan fftPlanR2C;
	ManagedArray<complex> temp; temp.init(nGdensePlanar);
	fftw_plan_with_nthreads(1); //Multiple simultaneous single threaded fourier transforms
	fftPlanR2C = fftw_plan_dft_r2c_2d(Sdense[jDir], Sdense[kDir], (double*)temp.data(), (fftw_complex*)temp.data(), FFTW_ESTIMATE);
	logPrintf("Done.\n");
	
	//Launch threads for initializing each plane perpendicular to axis:
	logPrintf("Computing truncated coulomb kernel ... "); logFlush();
	std::vector<CoulombKernelWire> ckwArr(nProcsAvailable);
	for(CoulombKernelWire& c: ckwArr)
	{	c.dense.init(nGdensePlanar*nArrays);
		c.Vc = data;
		c.Vc_RRT = data_RRT;
		c.fftPlanR2C = fftPlanR2C;
		//Copy geometry definitions:
		c.R = R; c.S = S; c.Sdense = Sdense;
		c.iDir = iDir; c.jDir = jDir; c.kDir = kDir;
		c.ws = &ws; c.sigma = sigma; c.omega = omega;
		c.arrayStride = nGdensePlanar; c.arrayStrideReal = 2*nGdensePlanar;
	}
	std::mutex mJobCount; int nPlanesDone = 0; //for job management
	threadLaunch(CoulombKernelWire::thread, 0, ckwArr.data(), 1+S[iDir]/2, &nPlanesDone, &mJobCount);
	fftw_destroy_plan(fftPlanR2C);
	logPrintf("Done.\n");
}
