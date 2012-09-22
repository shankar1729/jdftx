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

#ifndef JDFTX_CORE_LOOPMACROS_H
#define JDFTX_CORE_LOOPMACROS_H

//! One thread of a loop over real space (see applyFunc_r_sub for example)
#define THREAD_rLoop(code) \
	size_t i=iStart; \
	vector3<int> iv( \
		i / (S[2]*S[1]), \
		(i/S[2]) % S[1], \
		i % S[2] ); \
	while(i<iStop) \
	{	\
		code \
		\
		i++; if(i==iStop) break; \
		iv[2]++; \
		if(iv[2]==S[2]) \
		{	iv[2]=0; \
			iv[1]++; \
			if(iv[1]==S[1]) \
			{	iv[1] = 0; \
				iv[0]++; \
			} \
		} \
	}

//! The GPU equivalent of THREAD_rLoop
#define COMPUTE_rIndices \
	vector3<int> iv( \
		threadIdx.z + zBlock * blockDim.z, \
		kernelIndex(y), \
		kernelIndex(x) ); \
	if(iv[0]>=S[0] || iv[1]>=S[1] || iv[2]>=S[2]) return; \
	size_t i = iv[2] + S[2]*size_t(iv[1] + S[1]*iv[0]);


//! One thread of a loop over G-space (see L_sub in PHLO.cpp for example)
#define THREAD_fullGspaceLoop(code) \
	size_t i=iStart; \
	vector3<int> iG( i / (S[2]*S[1]), (i/S[2]) % S[1], i % S[2] ); \
	for(int j=0; j<3; j++) if(2*iG[j]>S[j]) iG[j]-=S[j]; \
	while(true) \
	{	\
		code \
		\
		i++; if(i==iStop) break; \
		iG[2]++; \
		if(2*iG[2]>S[2]) iG[2]-=S[2]; \
		if(iG[2]==0) \
		{	iG[1]++; \
			if(2*iG[1]>S[1]) iG[1]-=S[1]; \
			if(iG[1]==0) \
			{	iG[0]++; \
				if(2*iG[0]>S[0]) iG[0]-=S[0]; \
			} \
		} \
	}

//! The GPU equivalent of THREAD_fullGspaceLoop
//! NOTE: x and z are swapped in kernel indices since x is contiguous on GPU
#define COMPUTE_fullGindices \
	vector3<int> iG( \
		zBlock * blockDim.z + threadIdx.z, \
		kernelIndex(y), \
		kernelIndex(x) ); \
	if(iG[0]>=S[0] || iG[1]>=S[1] || iG[2]>=S[2]) return; \
	size_t i = iG[2] + S[2]*size_t(iG[1] + S[1]*iG[0]); \
	for(int j=0; j<3; j++) if(2*iG[j]>S[j]) iG[j]-=S[j];


//! One thread of a loop over symmetry-reduced G-space (see applyFuncGsq_sub in Operators.h for example)
#define THREAD_halfGspaceLoop(code) \
	int size2 = S[2]/2+1; \
	size_t i=iStart; \
	vector3<int> iG( i / (size2*S[1]), (i/size2) % S[1], i % size2 ); \
	for(int j=0; j<3; j++) if(2*iG[j]>S[j]) iG[j]-=S[j]; \
	while(i<iStop) \
	{	\
		code \
		\
		i++; if(i==iStop) break; \
		iG[2]++; \
		if(iG[2]==size2) \
		{	iG[2]=0; \
			iG[1]++; \
			if(2*iG[1]>S[1]) iG[1]-=S[1]; \
			if(iG[1]==0) \
			{	iG[0]++; \
				if(2*iG[0]>S[0]) iG[0]-=S[0]; \
			} \
		} \
	}

//! The GPU equivalent of THREAD_halfGspaceLoop
//! NOTE: x and z are swapped in kernel indices since x is contiguous on GPU
#define COMPUTE_halfGindices \
	int size2 = S[2]/2+1; \
	vector3<int> iG( \
		threadIdx.z + zBlock * blockDim.z, \
		kernelIndex(y), \
		kernelIndex(x) ); \
	if(iG[0]>=S[0] || iG[1]>=S[1] || iG[2]>=size2) return; \
	size_t i = iG[2] + size2*size_t(iG[1] + S[1]*iG[0]); \
	for(int j=0; j<3; j++) if(2*iG[j]>S[j]) iG[j]-=S[j];

//! Determine if this is a nyquist component
//! NOTE: no component is nyquist for odd sizes in this convention
#define IS_NYQUIST ( (!(2*iG[0]-S[0])) | (!(2*iG[1]-S[1])) | (!(2*iG[2]-S[2])) )

#endif // JDFTX_CORE_LOOPMACROS_H
