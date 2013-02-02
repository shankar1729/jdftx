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


#ifndef JDFTX_CORE_GPUKERNELUTILS_H
#define JDFTX_CORE_GPUKERNELUTILS_H

#include <algorithm>
#include <cuda_runtime.h>
#include <driver_types.h>
#include <vector_types.h>
#include <core/vector3.h>

//! @file GpuKernelUtils.h
//! @brief Common utility functions/macros for the gpu kernels and launchers in the .cu files

//! Base-class for launch configuration for gpu kernels
struct GpuLaunchConfig
{	cudaFuncAttributes attr; //!< attributes of the function
	cudaDeviceProp prop; //!< properties of the currnetly running device

	//! Initialize the device and function properties
	template<typename GpuKernel> GpuLaunchConfig(GpuKernel* gpuKernel)
	{	int iDevice; cudaGetDevice(&iDevice);
		cudaGetDeviceProperties(&prop, iDevice);
		cudaFuncGetAttributes(&attr, gpuKernel);
	}
};


//Get the logical index of the kernel (dir is x, y, or z)
#define kernelIndex(dir) (blockIdx.dir * blockDim.dir + threadIdx.dir)

//Get the logical 1D index, even if the grid is 2D (required for very large 1D kernels)
#define kernelIndex1D() ((blockIdx.y*gridDim.x+blockIdx.x) * blockDim.x + threadIdx.x)


//! 1D launch configuration
struct GpuLaunchConfig1D : public GpuLaunchConfig
{	dim3 nPerBlock; //!< dimension of block
	dim3 nBlocks; //!< dimension of grid (note nBlocks could be 3D for really large kernels)

	//! Set up blocks and grid for a 1D operation over N data points
	template<typename GpuKernel> GpuLaunchConfig1D(GpuKernel* gpuKernel, int N)
	: GpuLaunchConfig(gpuKernel),
	nPerBlock(attr.maxThreadsPerBlock,1,1),
	nBlocks(ceildiv(N, int(nPerBlock.x)),1,1)
	{	//If the grid is too big, make it 2D:
		while(int(nBlocks.x) > prop.maxGridSize[0])
		{	nBlocks.x = ceildiv(int(nBlocks.x),2);
			nBlocks.y *= 2;
		}
	}
};

//! 3D launch configuration
struct GpuLaunchConfig3D : public GpuLaunchConfig
{	dim3 nPerBlock; //!< dimension of block
	dim3 nBlocks; //!< dimension of grid (note nBlocks could be 3D for really large kernels)
	int zBlockMax; //!< Grids are 2D, so need to loop over last dim

	//! Set up blocks and grid for a 1D operation over N data points
	template<typename GpuKernel> GpuLaunchConfig3D(GpuKernel* gpuKernel, vector3<int> S)
	: GpuLaunchConfig(gpuKernel)
	{	// Try to minimize zBlockMax and maximize block size within constraint:
		zBlockMax = ceildiv(S[0], std::min(attr.maxThreadsPerBlock, prop.maxThreadsDim[2]));
		nPerBlock.z = ceildiv(S[0], zBlockMax);
		// For the chosen z configuration, maximize x block size within constraint
		int maxBlockXY = attr.maxThreadsPerBlock/nPerBlock.z;
		nBlocks.x = ceildiv(S[2], std::min(maxBlockXY,prop.maxThreadsDim[0]));
		nPerBlock.x = ceildiv(S[2], int(nBlocks.x));
		// For the chosen x and z configuration, maximize y block size within constraint
		int maxBlockY = attr.maxThreadsPerBlock/(nPerBlock.z*nPerBlock.x);
		nBlocks.y = ceildiv(S[1], std::min(maxBlockY,prop.maxThreadsDim[1]));
		nPerBlock.y = ceildiv(S[1], int(nBlocks.y));
	}
};

//! 3D launch configuration for symmetry-reduced G-space loops (z dimension folded for real data sets)
struct GpuLaunchConfigHalf3D : public GpuLaunchConfig3D
{	//!Just use the above after reducing the z-dimension to half
	template<typename GpuKernel> GpuLaunchConfigHalf3D(GpuKernel* gpuKernel, vector3<int> S)
	: GpuLaunchConfig3D(gpuKernel, vector3<int>(S[0], S[1], S[2]/2+1))
	{
	}
};

//! Check for gpu errors and print a useful message (implemented in GpuUtils.cpp)
void gpuErrorCheck();

#endif // JDFTX_CORE_GPUKERNELUTILS_H
