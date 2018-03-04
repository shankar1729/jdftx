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

#ifndef JDFTX_CORE_GPUUTIL_H
#define JDFTX_CORE_GPUUTIL_H

#ifdef GPU_ENABLED

#include <cstdio>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <vector>

//! @addtogroup Utilities
//! @{

//! @file GpuUtil.h

extern cublasHandle_t cublasHandle; //!< global handle to cublas (defined in GpuUtil.cpp)
#ifdef CUSOLVER_ENABLED
#include <cusolverDn.h>
extern cusolverDnHandle_t cusolverHandle;  //!< global handle to cusolverDn (defined in GpuUtil.cpp)
#endif

//! Must be called before any GPU use (preferably from main(), see #isGpuMine)
//! If mpiHostGpu (group of GPU processes on the same node) is specified, divide compatible GPUs amongst processes on same node, else select one with max memory
//! nGPUs returns the number of physical GPUs used (fraction if a GPU is shared with other processes)
//! Returns false on failure to find a suitable GPU
bool gpuInit(FILE* fpLog=stdout, const class MPIUtil* mpiHostGpu=0, double* nGPUs=0);

//! Only the thread that called gpuInit() is allowed to use GPU resources
//! This function will return true only on the one thread that rules the gpu.
bool isGpuMine();

//! Check for gpu errors, and if any, abort with a useful messsage
extern void gpuErrorCheck();

#endif //GPU_ENABLED


//! A utility function to eliminate ugly \#ifdef's when not required
inline bool isGpuEnabled()
{
	#ifdef GPU_ENABLED
	return true;
	#else
	return false;
	#endif
}

//! @}
#endif // JDFTX_CORE_GPUUTIL_H
