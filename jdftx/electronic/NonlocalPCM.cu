/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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
#include <core/LoopMacros.h>
#include <electronic/NonlocalPCM_internal.h>

__global__
void cavitationEnergy_kernel(int N, double* Earr,
	vector3<const double*> Dshape, symmMatrix3<const double*> DDshape,
	vector3<double*> grad_Dshape, symmMatrix3<double*> grad_DDshape)
{	int i = kernelIndex1D();
	if(i<N) Earr[i] = cavitationEnergy_calc(i, Dshape, DDshape, grad_Dshape, grad_DDshape);
}
void cavitationEnergy_gpu(int N, double* Earr,
	vector3<const double*> Dshape, symmMatrix3<const double*> DDshape,
	vector3<double*> grad_Dshape, symmMatrix3<double*> grad_DDshape)
{
	GpuLaunchConfig1D glc(cavitationEnergy_kernel, N);
	cavitationEnergy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, Earr, Dshape, DDshape, grad_Dshape, grad_DDshape);
}
