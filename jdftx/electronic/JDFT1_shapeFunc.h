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

#ifndef JDFTX_ELECTRONIC_JDFT1_SHAPEFUNC_H
#define JDFTX_ELECTRONIC_JDFT1_SHAPEFUNC_H

//The erfc 'shape function' of JDFT1 and its gradient at one point in space (shared between CPU and GPU)

#include <core/scalar.h>

__hostanddev__ void JDFT1_shapeFunc_sub(int i, const double* nCavity, double* shape, const double nc, const double sigma)
{	shape[i] = erfc(sqrt(0.5)*log(fabs(nCavity[i])/nc)/sigma)*0.5;
}

__hostanddev__ void JDFT1_shapeFunc_grad_sub(int i, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	grad_nCavity[i] = (-1.0/(nc*sigma*sqrt(2*M_PI))) * grad_shape[i]
		* exp(0.5*(pow(sigma,2) - pow(log(fabs(nCavity[i])/nc)/sigma + sigma, 2)));
}

//GPU kernel/launchers that loop over the above functions (implemented in JDFT1_shapeFunc.cu)
#ifdef GPU_ENABLED
void JDFT1_shapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma);
void JDFT1_shapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma);
#endif


#endif // JDFTX_ELECTRONIC_JDFT1_SHAPEFUNC_H
