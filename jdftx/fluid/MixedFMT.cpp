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

#include <fluid/MixedFMT.h>
#include <fluid/MixedFMT_internal.h>
#include <core/VectorField.h>

//Compute the tensor weighted density (threaded/gpu):
inline void tensorKernel_sub(size_t iStart, size_t iStop, vector3<int> S, const matrix3<> G,
	const complex* nTilde, tensor3<complex*> mTilde)
{	THREAD_halfGspaceLoop( tensorKernel_calc(i, iG, IS_NYQUIST, G, nTilde, mTilde); )
}
#ifdef GPU_ENABLED
void tensorKernel_gpu(vector3<int> S, const matrix3<> G,
	const complex* nTilde, tensor3<complex*> mTilde);
#endif
TensorFieldTilde tensorKernel(const ScalarFieldTilde& nTilde)
{	const GridInfo& gInfo = nTilde->gInfo;
	TensorFieldTilde mTilde(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	tensorKernel_gpu(gInfo.S, gInfo.G, nTilde->dataGpu(), mTilde.dataGpu());
	#else
	threadLaunch(tensorKernel_sub, gInfo.nG, gInfo.S, gInfo.G, nTilde->data(), mTilde.data());
	#endif
	return mTilde;
}

//Propagate gradient w.r.t tensor weighted density (threaded/gpu):
inline void tensorKernel_grad_sub(size_t iStart, size_t iStop, vector3<int> S, const matrix3<> G,
	tensor3<const complex*> grad_mTilde, complex* grad_nTilde)
{	THREAD_halfGspaceLoop( tensorKernel_grad_calc(i, iG, IS_NYQUIST, G, grad_mTilde, grad_nTilde); )
}
#ifdef GPU_ENABLED
void tensorKernel_grad_gpu(vector3<int> S, const matrix3<> G,
	tensor3<const complex*> grad_mTilde, complex* grad_nTilde);
#endif
ScalarFieldTilde tensorKernel_grad(const TensorFieldTilde& grad_mTilde)
{	const GridInfo& gInfo = grad_mTilde[0]->gInfo;
	ScalarFieldTilde grad_nTilde(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	tensorKernel_grad_gpu(gInfo.S, gInfo.G, grad_mTilde.const_dataGpu(), grad_nTilde->dataGpu());
	#else
	threadLaunch(tensorKernel_grad_sub, gInfo.nG, gInfo.S, gInfo.G, grad_mTilde.const_data(), grad_nTilde->data());
	#endif
	return grad_nTilde;
}


#ifdef GPU_ENABLED
void phiFMT_gpu(int N, double* phiArr,
	const double *n0arr, const double *n1arr, const double *n2arr, const double *n3arr,
	vector3<const double*> n1vArr, vector3<const double*> n2vArr, tensor3<const double*> n2mArr,
	double *grad_n0arr, double *grad_n1arr, double *grad_n2arr, double *grad_n3arr,
	vector3<double*> grad_n1vArr, vector3<double*> grad_n2vArr, tensor3<double*> grad_n2mArr);
#endif


double PhiFMT(const ScalarField& n0, const ScalarField& n1, const ScalarField& n2,
	const ScalarFieldTilde& n3tilde, const ScalarFieldTilde& n1vTilde, const ScalarFieldTilde& n2mTilde,
	ScalarField& grad_n0, ScalarField& grad_n1, ScalarField& grad_n2,
	ScalarFieldTilde& grad_n3tilde, ScalarFieldTilde& grad_n1vTilde, ScalarFieldTilde& grad_n2mTilde)
{
	const GridInfo& gInfo = n0->gInfo;
	ScalarField n3 = I(n3tilde);
	VectorField n1v = I(gradient(n1vTilde));
	VectorField n2v = I(gradient(-n3tilde));
	TensorField n2m = I(tensorKernel(n2mTilde));

	ScalarField grad_n3; VectorField grad_n1v, grad_n2v; TensorField grad_n2m;
	nullToZero(grad_n0, gInfo); nullToZero(grad_n1, gInfo); nullToZero(grad_n2, gInfo); nullToZero(grad_n3, gInfo);
	nullToZero(grad_n1v, gInfo); nullToZero(grad_n2v, gInfo); nullToZero(grad_n2m, gInfo);

	double result;
	#ifdef GPU_ENABLED
	{	ScalarField phiArr(ScalarFieldData::alloc(gInfo, true));
		phiFMT_gpu(gInfo.nr, phiArr->dataGpu(),
			n0->dataGpu(), n1->dataGpu(), n2->dataGpu(), n3->dataGpu(), n1v.const_dataGpu(), n2v.const_dataGpu(), n2m.const_dataGpu(),
			grad_n0->dataGpu(), grad_n1->dataGpu(), grad_n2->dataGpu(), grad_n3->dataGpu(),
			grad_n1v.dataGpu(), grad_n2v.dataGpu(), grad_n2m.dataGpu());
		result = gInfo.dV * sum(phiArr);
	}
	#else
	result = gInfo.dV*threadedAccumulate(phiFMT_calc, gInfo.nr,
			n0->data(), n1->data(), n2->data(), n3->data(), n1v.const_data(), n2v.const_data(), n2m.const_data(),
			grad_n0->data(), grad_n1->data(), grad_n2->data(), grad_n3->data(),
			grad_n1v.data(), grad_n2v.data(), grad_n2m.data());
	#endif
	n3=0; n1v=0; n2v=0; n2m=0; //no longer need these weighted densities (clean up)

	grad_n2mTilde += tensorKernel_grad(Idag(grad_n2m)); grad_n2m=0;
	grad_n1vTilde -= divergence(Idag(grad_n1v)); grad_n1v=0;
	grad_n3tilde += ( Idag(grad_n3) + divergence(Idag(grad_n2v)) ); grad_n3=0; grad_n2v=0;
	return result;
}

double phiFMTuniform(double n0, double n1, double n2, double n3,
	double& grad_n0, double& grad_n1, double& grad_n2, double& grad_n3)
{
	double zero=0.0, constzero=0.0;
	std::vector<double*> zeroArr(5,&zero);
	std::vector<const double*> constZeroArr(5,&constzero); //dummy arrays for zero vector and tensor weighted densities

	return phiFMT_calc(0, &n0, &n1, &n2, &n3, constZeroArr, constZeroArr, constZeroArr,
		&grad_n0, &grad_n1, &grad_n2, &grad_n3, zeroArr, zeroArr, zeroArr);
}

#ifdef GPU_ENABLED
void phiBond_gpu(int N, double Rhm, double scale, double* phiArr,
	const double *n0arr, const double *n2arr, const double *n3arr, vector3<const double*> n2vArr,
	double *grad_n0arr, double *grad_n2arr, double *grad_n3arr, vector3<double*> grad_n2vArr);
#endif

double PhiBond(double Rhm, double scale, const ScalarField& n0mol, const ScalarField& n2, const ScalarFieldTilde& n3tilde,
	ScalarField& grad_n0mol, ScalarField& grad_n2, ScalarFieldTilde& grad_n3tilde)
{
	const GridInfo& gInfo = n0mol->gInfo;
	//Compute n3 and n2v in real space from n3tilde:
	ScalarField n3 = I(n3tilde);
	VectorField n2v = I(gradient(-n3tilde));
	//Bonding correction and gradient:
	ScalarField grad_n3; VectorField grad_n2v;
	nullToZero(grad_n0mol, gInfo);
	nullToZero(grad_n2, gInfo);
	nullToZero(grad_n3, gInfo);
	nullToZero(grad_n2v, gInfo);
	#ifdef GPU_ENABLED
	ScalarField phiArr(ScalarFieldData::alloc(gInfo, true));
	phiBond_gpu(gInfo.nr, Rhm, scale, phiArr->dataGpu(),
		n0mol->dataGpu(), n2->dataGpu(), n3->dataGpu(), n2v.const_dataGpu(),
		grad_n0mol->dataGpu(), grad_n2->dataGpu(), grad_n3->dataGpu(), grad_n2v.dataGpu());
	double result = gInfo.dV * sum(phiArr);
	#else
	double result = gInfo.dV * threadedAccumulate(phiBond_calc, gInfo.nr, Rhm, scale,
			n0mol->data(), n2->data(), n3->data(), n2v.const_data(),
			grad_n0mol->data(), grad_n2->data(), grad_n3->data(), grad_n2v.data());
	#endif
	n3=0; n2v=0; //no longer need these weighted densities (clean up)
	//Propagate grad_n2v and grad_n3 to grad_n3tilde:
	grad_n3tilde += ( Idag(grad_n3) + divergence(Idag(grad_n2v)) );
	return result;
}

double phiBondUniform(double Rhm, double scale, double n0mol, double n2, double n3,
	double& grad_n0mol, double& grad_n2, double& grad_n3)
{	double zero=0.0, constzero=0.0;
	std::vector<double*> zeroArr(3,&zero);
	std::vector<const double*> constZeroArr(3,&constzero); //dummy arrays for zero vector weighted densities
	return phiBond_calc(0, Rhm, scale, &n0mol, &n2, &n3, constZeroArr, &grad_n0mol, &grad_n2, &grad_n3, zeroArr);
}
