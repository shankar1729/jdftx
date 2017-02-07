/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <fluid/TranslationOperator.h>
#include <core/Operators.h>
#include <core/Util.h>
#include <algorithm>

//! Special sparse matrix format for storing the linear spline translation (for a given displacement vector)
struct LsplineMatrix
{
	const GridInfo& gInfo;
	
	//Sparse Matrix representing S:
	int* index; //!< list of target indices for each source point (in range 0 to gInfo.S-2)
	double* weight; //!< right-side interpolation weight for each source point
	
	void S_axpy(double alpha, const ScalarField& x, ScalarField& y) const
	{	assert(x.gInfo == &gInfo);
		nullToZero(y, gInfo);
		const double* xData = x.data();
		double* yData = y.data();
		for(int j=0; j<gInfo.S; j++)
		{	int i = index[j];
			double w = weight[j];
			int iNext = w ? i+1 : i;
			yData[j] += alpha * ((1.-w) * xData[i] + w * xData[iNext]);
		}
	}
	
	void Sdag_axpy(double alpha, const ScalarField& x, ScalarField& y) const
	{	assert(x.gInfo == &gInfo);
		nullToZero(y, gInfo);
		const double* xData = x.data();
		double* yData = y.data();
		for(int j=0; j<gInfo.S; j++)
		{	int i = index[j];
			double w = weight[j];
			int iNext = w ? i+1 : i;
			double alpha_xj = alpha * xData[j];
			yData[i] += alpha_xj * (1.-w);
			yData[iNext] += alpha_xj * w;
		}
	}
	
	LsplineMatrix(const GridInfo& gInfo, const vector3<>& a) : gInfo(gInfo)
	{	index = new int[gInfo.S];
		weight = new double[gInfo.S];
		for(int j=0; j<gInfo.S; j++)
		{	double r = 0.;
			switch(gInfo.coord)
			{	case GridInfo::Spherical:   r = (a + vector3<>(0,0,gInfo.r[j])).length(); break;
				case GridInfo::Cylindrical: r = hypot(a[1], gInfo.r[j]+a[2]); break;
				case GridInfo::Planar:      r = gInfo.r[j]+a[2]; break;
			}
			//reduce r to within [0, rMax]
			while(r<0 || r>gInfo.rMax)
			{	if(r<0) r=-r;
				if(r>gInfo.rMax) r=2*gInfo.rMax-r;
			}
			//Find the interval to which r belongs, and assign index and weights appropriately:
			int i = 0;
			double w = 0.;
			if(gInfo.S > 1) //All data points to the 0th element in the special case of S=1
			{	if(r <= gInfo.r[0]) { i=0; w=0.; }
				else if(r >= gInfo.r[gInfo.S-1]) { i=gInfo.S-2; w=1.; }
				else
				{	i = int(std::upper_bound(gInfo.r.begin(), gInfo.r.end(), r) - gInfo.r.begin()) - 1;
					w = (r - gInfo.r[i]) / (gInfo.r[i+1] - gInfo.r[i]);
				}
			}
			index[j] = i;
			weight[j] = w;
		}
	}
	
	~LsplineMatrix()
	{	delete[] index;
		delete[] weight;
	}
};

TranslationOperatorLspline::TranslationOperatorLspline(const GridInfo& gInfo) : gInfo(gInfo)
{
}

void TranslationOperatorLspline::S_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const
{	getMatrix(a)->S_axpy(alpha, x, y);
}

void TranslationOperatorLspline::Sdag_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const
{	getMatrix(a)->Sdag_axpy(alpha, x, y);
}

std::shared_ptr<LsplineMatrix> TranslationOperatorLspline::getMatrix(const vector3<>& a) const
{	
	TranslationOperatorLspline* mThis = (TranslationOperatorLspline*)this; //modifiable copy of this
	
	mThis->lsplineMutex.lock();
	auto lsMatIter = lsplineMatrix.find(a);
	mThis->lsplineMutex.unlock();
	
	if(lsplineMatrix.find(a) == lsplineMatrix.end())
	{	//transltaion matrix for this displacement not yet computed
		auto lsMat = std::make_shared<LsplineMatrix>(gInfo, a);
		mThis->lsplineMutex.lock();
		mThis->lsplineMatrix[a] = lsMat;
		mThis->lsplineMutex.unlock();
		return lsMat;
	}
	else return lsMatIter->second; //previously computed
}
