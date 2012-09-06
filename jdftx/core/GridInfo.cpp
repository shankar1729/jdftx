/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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
#include <algorithm>

GridInfo::GridInfo():Gmax(0),initialized(false)
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
		fftw_destroy_plan(planForwardMulti);
		fftw_destroy_plan(planInverseMulti);
		fftw_destroy_plan(planForwardInPlaceMulti);
		fftw_destroy_plan(planInverseInPlaceMulti);
		fftw_destroy_plan(planCtoRmulti);
		fftw_destroy_plan(planRtoCmulti);
		fftw_cleanup_threads();
		fftw_cleanup();
		#endif
	}
}

void GridInfo::initialize()
{
	this->~GridInfo(); //cleanup previously initialized quantities

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

	logPrintf("\n---------- Initializing the Grid ----------\n");
	logPrintf("R = \n"); R.print(globalLog, "%10lg ");
	logPrintf("G =\n"); G.print(globalLog, "%10lg ");
	logPrintf("unit cell volume = %lg\n", detR);
	
	//Choose / verify validity of sample count S (fftbox size)
	vector3<int> Smin;
	if(Gmax) //Determine minimum sampling corresponding to Gmax:
	{	//The integer coefficient of the nyquist component (half S) in each direction
		//when multiplied by the spacing between the reciprocal lattice plane spacings
		//must just exceed twice Gmax for the density sphere to be inscribed:
		for(int k=0; k<3; k++) Smin[k] = std::max(1, 4*int(Gmax * R.column(k).length() / (2*M_PI)));
		logPrintf("Minimum fftbox size, Smin = "); Smin.print(globalLog, " %d ");
	}
	for(int k=0; k<3; k++)
	{	if(S[k]<=0) //this component not set, pick minimal FFT-suitable size
		{	assert(Smin[k]>0); //One of Gmax or S must be specified
			S[k] = Smin[k];
			while(!fftSuitable(S[k])) S[k]+=2; //move through even numbers
		}
		else //check validity:
		{	if(S[k]<Smin[k])
				die("Specified fftbox dimension S[%d] = %d < %d = Smin[%d] for the G-sphere.\n",
					k, S[k], Smin[k], k)
		}
	}
	logPrintf("Chosen fftbox size, S = "); S.print(globalLog, " %d ");

	//Sample size dependent quantities:
	nr = S[0] * S[1] * S[2];
	nG = S[0] * S[1] * (S[2]/2+1);
	dV = detR/nr;
	for(int k=0; k<3; k++) h[k] = R.column(k)/S[k];

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
	
	//Check for previously saved wisdom (this file should not be copied from one machine to another)
	int systemWisdom = fftw_import_system_wisdom();
	FILE* fp = fopen(".fftw-wisdom", "r");
	if(fp) { fftw_import_wisdom_from_file(fp); fclose(fp); }
	else if(!systemWisdom) logPrintf("No local or system FFTW wisdom found, planning might take a while ...\n");

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
