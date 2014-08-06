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

#include <core/Thread.h>

#include <float.h>
#include <string.h>
#include <stdio.h>

#if defined(MKL_PROVIDES_BLAS) || defined(MKL_PROVIDES_FFT)
#include <mkl.h>
#endif

int getPhysicalCores()
{	FILE* pp = popen("awk '$1==\"physical\" && $2==\"id\" && !ID[$4] { PROCS++; ID[$4]=1; } $1=\"cpu\" && $2==\"cores\" { CORES=$4; }  END { print PROCS*CORES }' /proc/cpuinfo", "r");
	if(pp)
	{	int nCores = 0;
		fscanf(pp, "%d", &nCores);
		fclose(pp);
		if(nCores) return nCores; //should work on Linux: physical cores ignoring hyperthreading
	}
	return sysconf(_SC_NPROCESSORS_ONLN); //POSIX compatible fallback (does not account for hyperthreading)
}

int nProcsAvailable = getPhysicalCores();
bool threadOperators = true;

bool shouldThreadOperators()
{	return threadOperators;
}

void suspendOperatorThreading()
{	threadOperators = false;
	#if defined(MKL_PROVIDES_BLAS) || defined(MKL_PROVIDES_FFT)
	mkl_set_num_threads(1);
	#endif
}

void resumeOperatorThreading()
{	threadOperators = true;
	#if defined(MKL_PROVIDES_BLAS) || defined(MKL_PROVIDES_FFT)
	mkl_set_num_threads(nProcsAvailable);
	#endif
}


AutoThreadCount::AutoThreadCount(int minThreads, int minStats, const char* name, FILE* fpLog)
: minThreads(minThreads),minStats(minStats),settled(false),name(0),fpLog(fpLog)
{	nMax = nProcsAvailable;
	if(minThreads>nMax) minThreads=nMax;
	time = new double[nMax];
	count = new int[nMax];
	for(int i=0; i<nMax; i++) { time[i]=DBL_MAX; count[i]=0; }
	if(name)
	{	this->name = new char[strlen(name)+1];
		strcpy(this->name, name);
	}
}

AutoThreadCount::~AutoThreadCount()
{	if(name)
	{	fprintf(fpLog, "\nTiming statistics for %s:\n", name);
		for(int i=0; i<nMax; i++)
			if(count[i])
				fprintf(fpLog, "\tAt %d threads: %.1lf us (%d tries)\n", i+1, time[i], count[i]);

		if(settled) fprintf(fpLog, "\tOptimum thread count was determined to be %d.\n", nOpt);
		else fprintf(fpLog, "\tNever settled on an optimum thread count.\n");

		delete[] name;
	}
	delete[] time;
	delete[] count;
}

int AutoThreadCount::getThreadCount()
{	if(!shouldThreadOperators()) return 1;
	if(settled) return nOpt; //optimum value has been found and previously stored
	for(int i=minThreads-1; i<nMax; i++)
	{	if(count[i]<minStats) return i+1; //not enough stats for this #threads yet, so get it!
		if(i && time[i]>time[i-1]) //threads becoming a burden
		{	nOpt = i; //= (i-1)+1
			settled = true; return nOpt;
		}
		//else: not found optimum yet, try higher thread count
	}
	//Searched through all thread counts, and found that higher thread count always seems to work better:
	nOpt = nMax;
	settled = true; return nOpt;
}

void AutoThreadCount::submitTime(int nThreads, double runTime)
{	int i = nThreads-1;
	time[i] = (time[i]*count[i] + runTime)/(count[i]+1);
	count[i]++;
}
