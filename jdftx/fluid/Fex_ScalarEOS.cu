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
#include <fluid/Fex_ScalarEOS_internal.h>

__global__
void evalJeffereyAustinEOS_kernel(int nr, const double* Nbar, double* Aex, double* Aex_N, double Vhs, const JeffereyAustinEOS_eval eval)
{	int i = kernelIndex1D();
	if(i<nr) eval(i, Nbar, Aex, Aex_N, Vhs);
}
void evalJeffereyAustinEOS_gpu(int nr, const double* Nbar, double* Aex, double* Aex_N, double Vhs, const JeffereyAustinEOS_eval& eval)
{	GpuLaunchConfig1D glc(evalJeffereyAustinEOS_kernel, nr);
	evalJeffereyAustinEOS_kernel<<<glc.nBlocks, glc.nPerBlock>>>(nr, Nbar, Aex, Aex_N, Vhs, eval);
}


__global__
void evalTaoMasonEOS_kernel(int nr, const double* Nbar, double* Aex, double* Aex_N, double Vhs, const TaoMasonEOS_eval eval)
{	int i = kernelIndex1D();
	if(i<nr) eval(i, Nbar, Aex, Aex_N, Vhs);
}
void evalTaoMasonEOS_gpu(int nr, const double* Nbar, double* Aex, double* Aex_N, double Vhs, const TaoMasonEOS_eval& eval)
{	GpuLaunchConfig1D glc(evalTaoMasonEOS_kernel, nr);
	evalTaoMasonEOS_kernel<<<glc.nBlocks, glc.nPerBlock>>>(nr, Nbar, Aex, Aex_N, Vhs, eval);
}

