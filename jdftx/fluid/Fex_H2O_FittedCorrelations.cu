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

#include <core/GpuKernelUtils.h>
#include <fluid/Fex_H2O_FittedCorrelations_internal.h>

__global__
void Fex_H20_FittedCorrelations_kernel(int nr, const double* NObar, const double* NHbar,
	double* Fex, double* Phi_NObar, double* Phi_NHbar)
{	int i = kernelIndex1D();
	if(i<nr)
		Fex[i] = Fex_H2O_FittedCorrelations_calc(i, NObar, NHbar, Phi_NObar, Phi_NHbar);
}
void Fex_H20_FittedCorrelations_gpu(int nr, const double* NObar, const double* NHbar,
	double* Fex, double* Phi_NObar, double* Phi_NHbar)
{	GpuLaunchConfig1D glc(Fex_H20_FittedCorrelations_kernel, nr);
	Fex_H20_FittedCorrelations_kernel<<<glc.nBlocks, glc.nPerBlock>>>(nr, NObar, NHbar, Fex, Phi_NObar, Phi_NHbar);
}

