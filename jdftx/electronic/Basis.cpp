/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#include <electronic/Everything.h>
#include <cstdio>
#include <cmath>

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#endif

Basis::Basis()
{	gInfo = 0;
	nbasis = 0;
	iGarr = NULL;
	index = NULL;
	#ifdef GPU_ENABLED
	iGarrGpu = NULL;
	indexGpu = NULL;
	#endif
	indexPref = NULL;
	ownsData = false;
}

Basis::~Basis()
{	if(ownsData)
	{	gInfo = 0;
		nbasis = 0;
		delete[] iGarr;
		delete[] index;
		#ifdef GPU_ENABLED
		cudaFree(iGarrGpu);
		cudaFree(indexGpu);
		#endif
		indexPref = NULL;
		ownsData = false;
	}
}

Basis::Basis(const Basis& basis)
{	*this = basis;
}

Basis& Basis::operator=(const Basis& basis)
{	gInfo = basis.gInfo;
	iInfo = basis.iInfo;
	nbasis = basis.nbasis;
	iGarr = basis.iGarr;
	index = basis.index;
	#ifdef GPU_ENABLED
	iGarrGpu = basis.iGarrGpu;
	indexGpu = basis.indexGpu;
	#endif
	indexPref = basis.indexPref;
	iGarrPref = basis.iGarrPref;
	head = basis.head;
	ownsData = false; //referenced all the other data, so don't own data
	return *this;
}


void Basis::setup(const GridInfo& gInfo, const IonInfo& iInfo, double Ecut, const vector3<> k)
{	//Find the indices within Ecut:
	vector3<int> iGbox; for(int i=0; i<3; i++) iGbox[i] = 1 + int(sqrt(2*Ecut) * gInfo.R.column(i).length() / (2*M_PI));
	std::vector< vector3<int> > iGvec;
	std::vector<int> indexVec;
	vector3<int> iG;
	for(iG[0]=-iGbox[0]; iG[0]<=iGbox[0]; iG[0]++)
		for(iG[1]=-iGbox[1]; iG[1]<=iGbox[1]; iG[1]++)
			for(iG[2]=-iGbox[2]; iG[2]<=iGbox[2]; iG[2]++)
				if(0.5*dot(iG+k, gInfo.GGT*(iG+k)) <= Ecut)
				{	iGvec.push_back(iG);
					indexVec.push_back(gInfo.fullGindex(iG));
				}
	setup(gInfo, iInfo, indexVec, iGvec);
	logPrintf("nbasis = %lu for k = ", nbasis); k.print(globalLog, " %6.3f ");
}

void Basis::setup(const GridInfo& gInfo, const IonInfo& iInfo, const std::vector<int>& indexVec)
{	//Compute the integer G-vectors for the specified indices:
	std::vector< vector3<int> > iGvec(indexVec.size());
	int stride1 = gInfo.S[2];
	int stride0 = gInfo.S[1] * stride1;
	for(unsigned j=0; j<indexVec.size(); j++)
	{	int index = indexVec[j];
		vector3<int> iG;
		iG[0] = index/stride0; index -= stride0*iG[0];
		iG[1] = index/stride1; index -= stride1*iG[1];
		iG[2] = index;
		for(int k=0; k<3; k++)
			if(2*iG[k] > gInfo.S[k])
				iG[k] -= gInfo.S[k];
		iGvec[j] = iG;
	}
	setup(gInfo, iInfo, indexVec, iGvec);
}



void Basis::setup(const GridInfo& gInfo, const IonInfo& iInfo,
	const std::vector<int>& indexVec, const std::vector< vector3<int> >& iGvec)
{
	this->gInfo = &gInfo;
	this->iInfo = &iInfo;
	
	nbasis = iGvec.size();
	iGarr = new vector3<int>[nbasis];
	index = new int[nbasis];
	memcpy(iGarr, &iGvec[0], sizeof(vector3<int>)*nbasis);
	memcpy(index, &indexVec[0], sizeof(int)*nbasis);

	#ifdef GPU_ENABLED
	//Copy index arrays over to GPU
	cudaMalloc(&indexGpu, sizeof(int)*nbasis);
	cudaMalloc(&iGarrGpu, sizeof(vector3<int>)*nbasis);
	gpuErrorCheck();
	cudaMemcpy(indexGpu, index, sizeof(int)*nbasis, cudaMemcpyHostToDevice);
	cudaMemcpy(iGarrGpu, iGarr, sizeof(vector3<int>)*nbasis, cudaMemcpyHostToDevice);
	gpuErrorCheck();
	indexPref = indexGpu;
	iGarrPref = iGarrGpu;
	#else
	indexPref = index;
	iGarrPref = iGarr;
	#endif

	//Initialize head:
	head.clear();
	for(size_t n=0; n<nbasis; n++)
		if(iGvec[n].length_squared() < 4) //selects 27 entries (basically [-1,+1]^3)
			head.push_back(n);
	
	ownsData = true;
}
