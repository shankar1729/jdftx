/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman, Deniz Gunceler

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

#include <electronic/PCM_internal.h>
#include <core/Operators.h>

//----------------------- The JDFT `shape function' and gradient ------------------

void pcmShapeFunc(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	threadedLoop(pcmShapeFunc_calc, N, nCavity, shape, nc, sigma);
}
void pcmShapeFunc_grad(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	threadedLoop(pcmShapeFunc_grad_calc, N, nCavity, grad_shape, grad_nCavity, nc, sigma);
}
#ifdef GPU_ENABLED
void pcmShapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma);
void pcmShapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma);
#endif
void pcmShapeFunc(const DataRptr& nCavity, DataRptr& shape, const double nc, const double sigma)
{	nullToZero(shape, nCavity->gInfo);
	callPref(pcmShapeFunc)(nCavity->gInfo.nr, nCavity->dataPref(), shape->dataPref(), nc, sigma);
}
void pcmShapeFunc_grad(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const double nc, const double sigma)
{	nullToZero(grad_nCavity, nCavity->gInfo);
	callPref(pcmShapeFunc_grad)(nCavity->gInfo.nr, nCavity->dataPref(), grad_shape->dataPref(), grad_nCavity->dataPref(), nc, sigma);
}
