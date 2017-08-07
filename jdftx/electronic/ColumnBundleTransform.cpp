/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <electronic/ColumnBundleTransform.h>
#include <electronic/ColumnBundle.h>
#include <electronic/Symmetries.h>
#include <core/LatticeUtils.h>
#include <core/BlasExtra.h>
#include <algorithm>

ColumnBundleTransform::BasisWrapper::BasisWrapper(const Basis& basis) : basis(basis)
{	//Determine bounds on iG:
	iGbox = vector3<int>();
	for(const vector3<int>& iG: basis.iGarr)
		for(int i=0; i<3; i++)
			iGbox[i] = std::max(iGbox[i], abs(iG[i]));
	//Initialize look-up table:
	pitch[2] = 1;
	pitch[1] = pitch[2] * (2*iGbox[2]+1);
	pitch[0] = pitch[1] * (2*iGbox[1]+1);
	table.assign(pitch[0] * (2*iGbox[0]+1), -1); //unused parts of box will contain -1
	size_t n=0;
	for(const vector3<int>& iG: basis.iGarr)
		table[dot(pitch,iG+iGbox)] = n++; //valid parts of box will contain corresponding index in basis
}

ColumnBundleTransform::ColumnBundleTransform(const vector3<>& kC, const Basis& basisC, const vector3<>& kD,
	const ColumnBundleTransform::BasisWrapper& basisDwrapper, int nSpinor, const SpaceGroupOp& sym, int invert, const matrix3<int>& super)
: basisC(basisC), basisD(basisDwrapper.basis), nSpinor(nSpinor), invert(invert)
{
	//Check k-point transformation and determine offset
	const matrix3<>& metricC = basisC.gInfo->RTR;
	assert(nrm2(metricC - (~sym.rot)*metricC*sym.rot) < symmThreshold * nrm2(metricC)); //check symmetry
	assert(abs(invert) == 1); //check inversion
	assert(nrm2(basisC.gInfo->R * super - basisD.gInfo->R) < symmThreshold * nrm2(basisD.gInfo->R)); //check supercell
	matrix3<int> affine = sym.rot * invert * super; //net affine transformation
	double offsetErr;
	vector3<int> offset = round(kC * affine - kD, &offsetErr);
	assert(offsetErr < symmThreshold);
	
	//Initialize index map:
	index.init(basisC.nbasis);
	int* indexPtr = index.data();
	for(const vector3<int>& iG_C: basisC.iGarr) //for each C recip lattice coords
	{	vector3<int> iG_D = iG_C * affine + offset; //corresponding D recip lattice coords
		*(indexPtr++) = basisDwrapper.table[dot(basisDwrapper.pitch, iG_D + basisDwrapper.iGbox)]; //use lookup table to get D index
	}
	assert(*std::min_element(index.begin(), index.end()) >= 0); //make sure all entries were found
	
	//Initialize translation phase (if necessary)
	if(sym.a.length_squared())
	{	phase.init(basisC.nbasis);
		complex* phasePtr = phase.data();
		for(const vector3<int>& iG_C: basisC.iGarr)
			*(phasePtr++) = cis((2*M_PI)*dot(iG_C + kC, sym.a));
	}
	
	//Initialize spinor transformation:
	switch(nSpinor)
	{	case 1: spinorRot = eye(1); break;
		case 2: 
		{	spinorRot = Symmetries::getSpinorRotation(~(basisC.gInfo->R * sym.rot * inv(basisC.gInfo->R)));
			if(invert<0)
			{	matrix sInvert = zeroes(2,2);
				sInvert.set(0,1, 1.);
				sInvert.set(1,0, -1.);
				spinorRot = conj(sInvert * spinorRot);
			}
			break;
		}
		default: assert(!"Invalid value for nSpinor");
	}
}

void ColumnBundleTransform::scatterAxpy(complex alpha, const ColumnBundle& C_C, int bC, ColumnBundle& C_D, int bD) const
{	//Check inputs:
	assert(C_C.colLength() == nSpinor*basisC.nbasis); assert(bC >= 0 && bC < C_C.nCols());
	assert(C_D.colLength() == nSpinor*basisD.nbasis); assert(bD >= 0 && bD < C_D.nCols());
	//Scatter:
	for(int sD=0; sD<nSpinor; sD++)
		for(int sC=0; sC<nSpinor; sC++)
			callPref(eblas_scatter_zaxpy)(index.nData(), alpha*spinorRot(sD,sC), index.dataPref(),
				C_C.dataPref() + C_C.index(bC, sC*C_C.basis->nbasis),
				C_D.dataPref() + C_D.index(bD, sD*C_D.basis->nbasis), invert<0,
				phase.dataPref(), false);
}

void ColumnBundleTransform::gatherAxpy(complex alpha, const ColumnBundle& C_D, int bD, ColumnBundle& C_C, int bC) const
{	//Check inputs:
	assert(C_C.colLength() == nSpinor*basisC.nbasis); assert(bC >= 0 && bC < C_C.nCols());
	assert(C_D.colLength() == nSpinor*basisD.nbasis); assert(bD >= 0 && bD < C_D.nCols());
	//Gather:
	matrix spinorRotInv = (invert<0) ? transpose(spinorRot) : dagger(spinorRot);
	for(int sD=0; sD<nSpinor; sD++)
		for(int sC=0; sC<nSpinor; sC++)
			callPref(eblas_gather_zaxpy)(index.nData(), alpha*spinorRotInv(sC,sD), index.dataPref(),
				C_D.dataPref() + C_D.index(bD, sD*C_D.basis->nbasis),
				C_C.dataPref() + C_C.index(bC, sC*C_C.basis->nbasis), invert<0,
				phase.dataPref(), !(invert < 0));
}

void ColumnBundleTransform::scatterAxpy(complex alpha, const ColumnBundle& C_C, ColumnBundle& C_D, int bDstart, int bDstep) const
{	for(int bC=0; bC<C_C.nCols(); bC++) scatterAxpy(alpha, C_C,bC, C_D,bDstart+bDstep*bC);
}

void ColumnBundleTransform::gatherAxpy(complex alpha, const ColumnBundle& C_D, int bDstart, int bDstep, ColumnBundle& C_C) const
{	for(int bC=0; bC<C_C.nCols(); bC++) gatherAxpy(alpha, C_D,bDstart+bDstep*bC, C_C,bC);
}

