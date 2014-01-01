/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler

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

#include <core/GridInfo.h>
#include <core/Util.h>
#include <core/Thread.h>
#include <core/Operators.h>
#include <core/LatticeUtils.h>
#include <algorithm>

#ifdef MKL_ENABLED
#include <fftw3_mkl.h>
#endif

const double GridInfo::maxAllowedStrain = 0.35;

GridInfo::GridInfo():Gmax(0),GmaxRho(0),initialized(false)
{
}

GridInfo::~GridInfo()
{
	if(initialized)
	{
		#ifdef GPU_ENABLED
		cufftDestroy(planZ2Z);
		cufftDestroy(planD2Z);
		cufftDestroy(planZ2D);
		cufftDestroy(planZ2Dcompat);
		#else
		fftw_destroy_plan(planForwardSingle);
		fftw_destroy_plan(planInverseSingle);
		fftw_destroy_plan(planForwardInPlaceSingle);
		fftw_destroy_plan(planInverseInPlaceSingle);
		fftw_destroy_plan(planCtoRsingle);
		fftw_destroy_plan(planRtoCsingle);
		if(nProcsAvailable > 1)
		{	fftw_destroy_plan(planForwardMulti);
			fftw_destroy_plan(planInverseMulti);
			fftw_destroy_plan(planForwardInPlaceMulti);
			fftw_destroy_plan(planInverseInPlaceMulti);
			fftw_destroy_plan(planCtoRmulti);
			fftw_destroy_plan(planRtoCmulti);
		}
		#endif
	}
}

void GridInfo::update()
{
	detR = fabs(det(R));
	
	RT = ~R;
	RTR = RT*R;
	invR = inv(R);
	invRT = inv(RT);
	invRTR = inv(RTR);

	G = (2.0*M_PI)*inv(R);
	GT = ~G;
	GGT = G*GT;
	invGGT = inv(GGT);
	
	if(nr) updateSdependent();
}

void GridInfo::updateSdependent()
{
	dV = detR/nr;
	for(int k=0; k<3; k++) h[k] = R.column(k)/S[k];
	
	//Recommended extent for radial functions:
	dGradial = 0.02;
	GmaxGrid = 0.0;
	{	vector3<int> c;
		for(c[0]=-1; c[0]<=1; c[0]+=2) for(c[1]=-1; c[1]<=1; c[1]+=2) for(c[2]=-1; c[2]<=1; c[2]+=2)
		{	vector3<> f; for(int k=0; k<3; k++) f[k] = c[k]*(S[k]/2);
			double g = sqrt(GGT.metric_length_squared(f));
			if(g>GmaxGrid) GmaxGrid=g;
		}
	}
	GmaxSphere = Gmax * (1.+maxAllowedStrain);
}

void GridInfo::printLattice()
{	logPrintf("R = \n"); R.print(globalLog, "%10lg ");
	logPrintf("unit cell volume = %lg\n", detR);
}

void GridInfo::printReciprocalLattice()
{	logPrintf("G =\n"); G.print(globalLog, "%10lg ");
}


void GridInfo::setLatticeVectors()
{	//Check manually specified lattice vectors
	if(latticeType == Manual)
	{	if(fabs(det(R)) < symmThresholdSq) throw string("Specified lattice vectors seem to be linearly dependent");
		return;
	}
	//Check parameter ranges:
	if(a <= 0.) throw string("Lattice dimension a must be positive.");
	if(b <= 0.) throw string("Lattice dimension b must be positive.");
	if(c <= 0.) throw string("Lattice dimension c must be positive.");
	if(alpha <= 0.) throw string("Lattice angle alpha must be positive.");
	if(beta  <= 0.) throw string("Lattice angle beta must be positive.");
	if(gamma <= 0.) throw string("Lattice angle gamma must be positive.");
	if(alpha >= 180.) throw string("Lattice angle alpha must be < 180 degrees.");
	if(beta  >= 180.) throw string("Lattice angle beta must be < 180 degrees.");
	if(gamma >= 180.) throw string("Lattice angle gamma must be < 180 degrees.");
	if(alpha > beta + gamma) throw string("Lattice angles must satisfy beta + gamma > alpha");
	if(beta > gamma + alpha) throw string("Lattice angles must satisfy gamma + alpha > beta");
	if(gamma > alpha + beta) throw string("Lattice angles must satisfy alpha + beta > gamma");
	//Compute base lattice vectors:
	double cosAlpha = cos(alpha * M_PI/180.);
	double cosBeta = cos(beta * M_PI/180.);
	double sinGamma, cosGamma; sincos(gamma * M_PI/180., &sinGamma, &cosGamma);
	R.set_col(0, a * vector3<>(1., 0., 0.));
	R.set_col(1, b * vector3<>(cosGamma, sinGamma, 0));
	R.set_col(2, c * vector3<>(cosBeta, (cosAlpha-cosBeta*cosGamma)/sinGamma, 0.));
	R(2,2) = sqrt(c*c - R(0,2)*R(0,2) - R(1,2)*R(1,2));
	//Set lattice vectors with appropriate modification:
	switch(latticeModification)
	{	case Simple: //No modification
			break;
		case BodyCentered:
		{	switch(latticeType)
			{	case Orthorhombic: case Tetragonal: case Cubic: break;
				default: throw string("Body-centered modification is valid only for Orthorhombic, Tetragonal or Cubic lattices");
			}
			R = R * 0.5 * matrix3<>
				(	-1, 1, 1,
					1, -1, 1,
					1, 1, -1 );
			break;
		}
		case BaseCentered:
		{	switch(latticeType)
			{	case Monoclinic: case Orthorhombic: break;
				default: throw string("Base-centered modification is valid only for Monoclinic or Orthorhombic lattices");
			}
			R = R * 0.5 * matrix3<>
				(	1, -1, 0,
					1, 1, 0,
					0, 0, 2 );
			break;
		}
		case FaceCentered:
		{	switch(latticeType)
			{	case Orthorhombic: case Cubic: break;
				default: throw string("Face-centered modification is valid only for Orthorhombic or Cubic lattices");
			}
			R = R * 0.5 * matrix3<>
				(	0, 1, 1,
					1, 0, 1,
					1, 1, 0 );
			break;
		}
	}
	for(int j=0; j<3; j++) for(int k=0; k<3; k++) if(fabs(R(j,k))<symmThresholdSq) R(j,k) = 0.; //Prevent annoying scientific notation printouts of almost zeros
}

//recursive function for symmetry-constrained S determination (see below)
void processSb(vector3<int>& Sb, int j, const matrix3<int>& ratios, vector3<bool>& dimsDone)
{	dimsDone[j] = true;
	for(int k=0; k<3; k++)
		if(j!=k && (ratios(j,k) || ratios(k,j)))
		{	if(Sb[k]) //verify consistency of pre-existing basis entry
			{	if( (ratios(j,k)*Sb[j])%Sb[k] || (ratios(k,j)*Sb[k])%Sb[j] ) //basis violates constraints between j and k
				{	logPrintf("WARNING: Could not find anisotropic solution for S commensurate with symmetries. Falling back to an isotropic solution.\n");
					Sb = vector3<int>(1,1,1);
				}
			}
			else //add k to this basis entry
			{	if(ratios(k,j))
				{	int x = Sb[j];
					Sb *= ratios(k,j); 
					Sb[k] = x * std::max(1, ratios(j,k));
				}
				else //ratios(j,k) is non-zero
				{	Sb[k] = ratios(j,k) * Sb[j];
				}
				processSb(Sb, k, ratios, dimsDone);
			}
		}
}

void GridInfo::initialize(bool skipHeader, const std::vector< matrix3<int> > sym)
{
	this->~GridInfo(); //cleanup previously initialized quantities
	
	update();
	if(!skipHeader) logPrintf("\n---------- Initializing the Grid ----------\n");
	printLattice();
	printReciprocalLattice();
	
	//Choose / verify validity of sample count S (fftbox size)
	vector3<int> Smin;
	bool autoS = (S[0]<=0) || (S[1]<=0) || (S[2]<=0); //whether to auto-compute sample count
	if(GmaxRho && GmaxRho < 2*Gmax) die("Charge density G-cutoff must be atleast twice that for wavefunctions.\n");
	if(Gmax || GmaxRho) //Determine minimum sampling corresponding to Gmax:
	{	//The integer coefficient of the nyquist component (half S) in each direction
		//when multiplied by the spacing between the reciprocal lattice plane spacings
		//must just exceed twice Gmax for the density sphere to be inscribed:
		for(int k=0; k<3; k++) Smin[k] = std::max(1, 4*int(std::max(Gmax,0.5*GmaxRho) * R.column(k).length() / (2*M_PI)));
		logPrintf("Minimum fftbox size, Smin = "); Smin.print(globalLog, " %d ");
	}
	if(autoS) //pick minimal FFT-suitable size
	{	//Determine constraints on S due to symmetries:
		matrix3<int> ratios;
		for(const matrix3<int>& m: sym)
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
					if(m(j,k))
						ratios(j,k) = gcd(ratios(j,k), abs(m(j,k)));
		//Construct integer basis of S's that satisfy these constraints:
		S = vector3<int>(0,0,0);
		vector3<bool> dimsDone(false,false,false); //dimensions yet to be covered by Sbasis
		for(int j=0; j<3; j++) if(!dimsDone[j])
		{	vector3<int> Sb; Sb[j] = 1;
			processSb(Sb, j, ratios, dimsDone);
			Sb = gcdReduce(Sb);
			//Check FFT suitability of the integers in Sb:
			bool isSuitable = true;
			for(int k=0; k<3; k++) if(Sb[k]) isSuitable = isSuitable && fftSuitable(Sb[k]);
			if(!isSuitable)
			{	logPrintf("WARNING: Symmetries require anisotropic S to include FFT-unsupported factors. Falling back to isotropic solution.\n");
				for(int k=0; k<3; k++) if(Sb[k]) Sb[k] = 1;
			}
			//For each basis entry, determine smallest fft-suitable scale factor that satisfies Smin constraint:
			int scaleSb = 0;
			for(int k=0; k<3; k++) if(Sb[k])
			{	int s = 2*ceildiv(Smin[k], 2*Sb[k]); //ensure even
				if(s > scaleSb) scaleSb = s;
			}
			while(!fftSuitable(scaleSb)) scaleSb += 2; //move through even numbers
			Sb *= scaleSb;
			S += Sb;
		}
	}
	else //Manually-specified sample count, only check validity:
	{	for(int k=0; k<3; k++)
		{	if(S[k]<Smin[k])
				die("Specified fftbox dimension S[%d] = %d < %d = Smin[%d] for the G-sphere.\n",
					k, S[k], Smin[k], k)
		}
	}
	logPrintf("Chosen fftbox size, S = "); S.print(globalLog, " %d ");

	//Sample size dependent quantities:
	nr = S[0] * S[1] * S[2];
	nG = S[0] * S[1] * (S[2]/2+1);
	updateSdependent();
	
	//Process division recommendations:
	irStart = (nr * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	irStop = (nr * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	iGstart = (nG * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	iGstop = (nG * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	
	logPrintf("Planning FFTs (this might take a while for a new big problem size) ... "); logFlush();

	#ifdef GPU_ENABLED
	//########## CUFFT Initialization ##############
	cufftPlan3d(&planZ2Z, S[0], S[1], S[2], CUFFT_Z2Z);
	cufftPlan3d(&planD2Z, S[0], S[1], S[2], CUFFT_D2Z);
	cufftPlan3d(&planZ2D, S[0], S[1], S[2], CUFFT_Z2D);
	cufftPlan3d(&planZ2Dcompat, S[0], S[1], S[2], CUFFT_Z2D);
	cufftSetCompatibilityMode(planZ2Dcompat, CUFFT_COMPATIBILITY_FFTW_ALL);
	gpuErrorCheck();
	
	#else
	//########## FFTW Initialization ##############
	#define PLANNER_FLAGS FFTW_MEASURE
	fftw_init_threads();
	#ifdef MKL_ENABLED
	fftw3_mkl.number_of_user_threads = nProcsAvailable; //Single-threaded transforms may be called from multiple threads
	#endif

	//Check for previously saved wisdom (this file should not be copied from one machine to another)
	int systemWisdom = fftw_import_system_wisdom();
	FILE* fp = fopen(".fftw-wisdom", "r");
	if(fp) { fftw_import_wisdom_from_file(fp); fclose(fp); }
	else if(!systemWisdom) logPrintf("\nNo local or system FFTW wisdom found, planning might take a while ... ");

	//Temp data for planning:
	fftw_complex* testData = (fftw_complex*)fftw_malloc(nr*sizeof(complex));
	fftw_complex* testData2 = (fftw_complex*)fftw_malloc(nr*sizeof(complex));
	
	//Single-threaded plans
	fftw_plan_with_nthreads(1);
	planInverseSingle = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData2, FFTW_BACKWARD, FFTW_MEASURE);
	planForwardSingle = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData2, FFTW_FORWARD, FFTW_MEASURE);
	planInverseInPlaceSingle = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData, FFTW_BACKWARD, FFTW_MEASURE);
	planForwardInPlaceSingle = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData, FFTW_FORWARD, FFTW_MEASURE);
	planRtoCsingle = fftw_plan_dft_r2c_3d(S[0], S[1], S[2], (double*)testData, testData2, PLANNER_FLAGS);
	planCtoRsingle = fftw_plan_dft_c2r_3d(S[0], S[1], S[2], testData, (double*)testData2, PLANNER_FLAGS);
	if(!planInverseSingle || !planForwardSingle
		|| !planInverseInPlaceSingle || !planForwardInPlaceSingle
		|| !planRtoCsingle || !planCtoRsingle)
		die("\nSingle-threaded FFTW planning failed in GridInfo::initialize()\n\n")
	
	//Multi-threaded plans:
	if(nProcsAvailable > 1)
	{	fftw_plan_with_nthreads(nProcsAvailable);
		planInverseMulti = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData2, FFTW_BACKWARD, FFTW_MEASURE);
		planForwardMulti = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData2, FFTW_FORWARD, FFTW_MEASURE);
		planInverseInPlaceMulti = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData, FFTW_BACKWARD, FFTW_MEASURE);
		planForwardInPlaceMulti = fftw_plan_dft_3d(S[0], S[1], S[2], testData, testData, FFTW_FORWARD, FFTW_MEASURE);
		planRtoCmulti = fftw_plan_dft_r2c_3d(S[0], S[1], S[2], (double*)testData, testData2, PLANNER_FLAGS);
		planCtoRmulti = fftw_plan_dft_c2r_3d(S[0], S[1], S[2], testData, (double*)testData2, PLANNER_FLAGS);
		if(!planInverseMulti || !planForwardMulti
			|| !planInverseInPlaceMulti || !planForwardInPlaceMulti
			|| !planRtoCmulti || !planCtoRmulti)
			die("\nMulti-threaded FFTW planning failed in GridInfo::initialize()\n\n")
	}
	else
	{	//Single-proc execution: copy over single threaded plans
		planInverseMulti = planInverseSingle;
		planForwardMulti = planForwardSingle;
		planInverseInPlaceMulti = planInverseInPlaceSingle;
		planForwardInPlaceMulti = planForwardInPlaceSingle;
		planRtoCmulti = planRtoCsingle;
		planCtoRmulti = planCtoRsingle;
	}
	fftw_free(testData);
	fftw_free(testData2);

	//Save wisdom in run directory for future use (copy .fftw-wisdom to /etc/fftw/wisdom for system-wide use)
	fp = fopen(".fftw-wisdom", "w");
	if(fp)
	{	fftw_export_wisdom_to_file(fp);
		fclose(fp);
	} //ignore if running in a read-only directory
	#endif

	logPrintf("Done.\n"); logFlush();
	initialized = true;
}
