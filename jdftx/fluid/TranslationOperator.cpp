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

#include <fluid/TranslationOperator.h>
#include <fluid/TranslationOperator_internal.h>
#include <core/Operators.h>


TranslationOperator::TranslationOperator(const GridInfo& gInfo) : gInfo(gInfo)
{
}

TranslationOperatorSpline::TranslationOperatorSpline(const GridInfo& gInfo, SplineType splineType)
: TranslationOperator(gInfo), splineType(splineType)
{
}

void constantSplineTaxpy_sub(size_t iStart, size_t iStop, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint)
{	THREAD_rLoop(constantSplineTaxpy_calc(i, iv, S, alpha, x, y, Tint);)
}
void linearSplineTaxpy_sub(size_t iStart, size_t iStop, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint, const vector3<> Tfrac)
{	THREAD_rLoop(linearSplineTaxpy_calc(i, iv, S, alpha, x, y, Tint, Tfrac);)
}
#ifdef GPU_ENABLED
void constantSplineTaxpy_gpu(const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint);
void linearSplineTaxpy_gpu(const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint, const vector3<> Tfrac);
#endif
void TranslationOperatorSpline::taxpy(const vector3<>& t, double alpha, const ScalarField& x, ScalarField& y) const
{	//Perform a gather with the inverse translation (hence negate t),
	//instead of scatter which is less efficient to parallelize
	vector3<> Tfrac = Diag(gInfo.S) * inv(gInfo.R) * (-t); //now in grid point units
	vector3<int> Tint;
	//Prepare output:
	nullToZero(y, gInfo);
	switch(splineType)
	{	case Constant:
		{	for(int k=0; k<3; k++)
			{	//round to nearest integer (and ensure symmetric rounding direction for transpose correctness):
				Tint[k] = int(copysign(floor(fabs(Tfrac[k])+0.5), Tfrac[k]));
				//reduce to positive first unit cell:
				Tint[k] = Tint[k] % gInfo.S[k];
				if(Tint[k]<0) Tint[k] += gInfo.S[k];
			}
			//Launch threads/gpu kernels:
			#ifdef GPU_ENABLED
			constantSplineTaxpy_gpu(gInfo.S, alpha*x->scale, x->dataGpu(false), y->dataGpu(), Tint);
			#else
			threadLaunch(constantSplineTaxpy_sub, gInfo.nr, gInfo.S, alpha*x->scale, x->data(false), y->data(), Tint);
			#endif
			break;
		}
		case Linear:
		{	for(int k=0; k<3; k++)
			{	//reduce to positive first unit cell:
				Tfrac[k] = fmod(Tfrac[k], gInfo.S[k]);
				if(Tfrac[k]<0) Tfrac[k] += gInfo.S[k];
				//separate integral and fractional parts:
				Tint[k] = int(floor(Tfrac[k]));
				Tfrac[k] -= Tint[k];
				Tint[k] = Tint[k] % gInfo.S[k];
			}
			//Launch threads/gpu kernels:
			#ifdef GPU_ENABLED
			linearSplineTaxpy_gpu(gInfo.S, alpha*x->scale, x->dataGpu(false), y->dataGpu(), Tint, Tfrac);
			#else
			threadLaunch(linearSplineTaxpy_sub, gInfo.nr, gInfo.S, alpha*x->scale, x->data(false), y->data(), Tint, Tfrac);
			#endif
			break;
		}
	}
}

TranslationOperatorFourier::TranslationOperatorFourier(const GridInfo& gInfo)
: TranslationOperator(gInfo)
{
}
inline void fourierTranslate_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> Gt, complex* xTilde)
{	THREAD_halfGspaceLoop( fourierTranslate_calc(i, iG, S, Gt, xTilde); )
}
#ifdef GPU_ENABLED //implemented in TranslationOperator.cu
void fourierTranslate_gpu(const vector3<int> S, const vector3<> Gt, complex* xTilde);
#endif
void TranslationOperatorFourier::taxpy(const vector3<>& t, double alpha, const ScalarField& x, ScalarField& y) const
{	ScalarFieldTilde xTilde = J(x);
	#ifdef GPU_ENABLED
	fourierTranslate_gpu(gInfo.S, gInfo.G*t, xTilde->dataGpu(false));
	#else
	threadLaunch(fourierTranslate_sub, gInfo.nG, gInfo.S, gInfo.G*t, xTilde->data(false));
	#endif
	y += alpha*I(xTilde, true);
}
