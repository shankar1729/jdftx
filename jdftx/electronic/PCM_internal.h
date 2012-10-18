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

#ifndef JDFTX_ELECTRONIC_PCM_INTERNAL_H
#define JDFTX_ELECTRONIC_PCM_INTERNAL_H

#include <core/scalar.h>

//----------- Common PCM functions (top level interface not seen by .cu files) ------------
#ifndef __in_a_cu_file__

#include <core/Data.h>

//! Compute the shape function (0 to 1) given the cavity-determining electron density
void pcmShapeFunc(const DataRptr& nCavity, DataRptr& shape, const double nc, const double sigma);

//!Compute derivative with respect to cavity-determining electron density, given derivative with respect to shape function
void pcmShapeFunc_grad(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const double nc, const double sigma);

#endif


//--------- Compute kernels (shared by CPU and GPU implementations) --------

//Cavity shape function and gradient
__hostanddev__ void pcmShapeFunc_calc(int i, const double* nCavity, double* shape, const double nc, const double sigma)
{	shape[i] = erfc(sqrt(0.5)*log(fabs(nCavity[i])/nc)/sigma)*0.5;
}
__hostanddev__ void pcmShapeFunc_grad_calc(int i, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	grad_nCavity[i] = (-1.0/(nc*sigma*sqrt(2*M_PI))) * grad_shape[i]
		* exp(0.5*(pow(sigma,2) - pow(log(fabs(nCavity[i])/nc)/sigma + sigma, 2)));
}


#endif // JDFTX_ELECTRONIC_PCM_INTERNAL_H
