/*-------------------------------------------------------------------
Copyright 2018 Ravishankar Sundararaman

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
{	FILE* pp = popen("awk '$1==\"physical\" && $2==\"id\" && !ID[$4] { PROCS++; ID[$4]=1; } $1=\"cpu\" && $2==\"cores\" { CORES=$4; }  END { print PROCS*CORES }' /proc/cpuinfo 2>/dev/null", "r");
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
	mkl_free_buffers();
	#ifndef THREADED_BLAS
	mkl_domain_set_num_threads(1, MKL_DOMAIN_BLAS); //Force single-threaded BLAS
	#endif
	#endif
}
