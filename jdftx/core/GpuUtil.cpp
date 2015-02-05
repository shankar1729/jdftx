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

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#include <core/Util.h>
#include <pthread.h>
#include <utility>
#include <cstdlib>

pthread_key_t gpuOwnerKey; //thread-local storage to identify thread that owns gpu
//NOTE: At the time of writing, c++0x threads implemented in g++, but not thread-local storage
//Using pthreads mechanism here, assuming that pthreads underly the c++0x threads
//This may not be true on Windows or for non-gcc compilers!

bool gpuInit(FILE* fpLog)
{	//Thread local storage to identify GPU owner thread
	pthread_key_create(&gpuOwnerKey, 0);
	pthread_setspecific(gpuOwnerKey, (const void*)1); //this will show up as 1 only on current thread

	//Selects the compatible GPU with maximum memory
	int nDevices, selectedDevice=-1; unsigned maxGlobalMem=0;
	cudaGetDeviceCount(&nDevices);
	for(int device=0; device<nDevices; device++)
	{	cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		std::pair<int,int> computeCap(prop.major, prop.minor);
		if(computeCap != std::make_pair(9999,9999) // not the emulation device
			&& computeCap >= std::make_pair(1,3) //compute capability >= 1.3 for double precision
			&& !prop.integrated) //reject on-board devices
		{
			fprintf(fpLog, "gpuInit: Found compatible cuda device %d '%s'\n", device, prop.name);
			if(prop.totalGlobalMem > maxGlobalMem)
			{	maxGlobalMem = prop.totalGlobalMem;
				selectedDevice = device;
			}
		}
	}
	if(selectedDevice < 0)
	{	fprintf(fpLog, "gpuInit: No compatible devices (>=1.3 compute capability, not on-board) found\n");
		return false;
	}
	fprintf(fpLog, "gpuInit: Selected device %d\n", selectedDevice);
	cudaSetDevice(selectedDevice);
	return true;
}

bool isGpuMine()
{	return bool(pthread_getspecific(gpuOwnerKey));
}

void gpuErrorCheck()
{	//cudaThreadSynchronize(); //NOTE: Uncomment this when trying to debug GPU kernel launches
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
	{	fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err));
		stackTraceExit(1);
	}
}


#endif //GPU_ENABLED
