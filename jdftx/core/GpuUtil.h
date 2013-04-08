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
#include <cublas.h>
#include <cufft.h>

//! @file GpuUtil.h

//! Must be called before any GPU use (preferably from main(), see #isGpuMine)
//! Returns false on failure to find a suitable GPU
bool gpuInit(FILE* fpLog=stdout);

//! Only the thread that called gpuInit() is allowed to use GPU resources
//! This function will return true only on the one thread that rules the gpu.
bool isGpuMine();

//! Check for gpu errors, and if any, abort with a useful messsage
extern void gpuErrorCheck();

#endif //GPU_ENABLED


//A utility function to eliminate ugly #ifdef's when not required
inline bool isGpuEnabled()
{
	#ifdef GPU_ENABLED
	return true;
	#else
	return false;
	#endif
}
#endif // JDFTX_CORE_GPUUTIL_H
