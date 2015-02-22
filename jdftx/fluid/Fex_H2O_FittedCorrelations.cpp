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

#include <fluid/Fex_H2O_FittedCorrelations_internal.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <electronic/operators.h>
#include <core/Units.h>
#include <core/Operators.h>

inline double COO_calc(double G)
{	static const double COO_A[6] = {-0.0271153, -0.0795576, 0.096648, -0.0291517, 0.0227052, -0.0109078};
	static const double COO_B[6] = {1.0569, 1.63853, 1.87815, 2.46787, 3.10047, 3.7162};
	static const double COO_C[6] = {0.0204818, 0.129513, 0.0736776, 0.0563636, 0.0844132, 0.0452023};
	double COO = 0.0;
	for(int j=0; j<6; j++) COO += COO_A[j] * exp(-pow(G-COO_B[j], 2) / COO_C[j]);
	return COO;
}

inline double COH_calc(double G)
{	static const double COH_A[6] = {-0.00801259, 0.0379836, -0.0380737, 0.0254576, -0.00291329, 0.00109967};
	static const double COH_B[6] = {1.19859, 1.72897, 2.15839, 2.72112, 3.29894, 3.7659};
	static const double COH_C[6] = {0.028028, 0.117534, 0.213336, 0.154157, 0.0265282, 0.0315518};
	double COH = 0.0;
	for(int j=0; j<6; j++) COH += COH_A[j] * exp(-pow(G-COH_B[j], 2) / COH_C[j]);
	return COH;
}

inline double CHH_calc(double G)
{	static const double CHH_A[2] = {-0.013959, 0.0295776};
	static const double CHH_B[2] = {1.88697, 2.53164};
	static const double CHH_C[2] = {0.104186, 0.0869848};
	double CHH = 0.0;
	for(int j=0; j<2; j++) CHH += CHH_A[j] * exp(-pow(G-CHH_B[j], 2) / CHH_C[j]);
	return CHH;
}

inline double fex_gauss_calc(double G)
{	const double r0 = 4.2027; //width of excess function correlations
	return exp(-pow(0.5*r0*G, 2));
}


Fex_H2O_FittedCorrelations::Fex_H2O_FittedCorrelations(const FluidMixture* fluidMixture, const FluidComponent* comp)
: Fex(fluidMixture,comp)
{
	if(fabs(T/Kelvin-298)>1)
		die("The FittedCorrelations functional is only valid at T=298K.\n")
	
	//Initialize the kernels:
	COO.init(0, gInfo.dGradial, gInfo.GmaxGrid, COO_calc);
	COH.init(0, gInfo.dGradial, gInfo.GmaxGrid, COH_calc);
	CHH.init(0, gInfo.dGradial, gInfo.GmaxGrid, CHH_calc);
	fex_gauss.init(0, gInfo.dGradial, gInfo.GmaxGrid, fex_gauss_calc);

	Citations::add("Fitted-Correlations water functional",
		"J. Lischner and T.A. Arias, J Phys Chem B. 114, 1946 (2010)");
}
Fex_H2O_FittedCorrelations::~Fex_H2O_FittedCorrelations()
{	COO.free();
	COH.free();
	CHH.free();
	fex_gauss.free();
}

#ifdef GPU_ENABLED
void Fex_H20_FittedCorrelations_gpu(int nr, const double* NObar, const double* NHbar,
	double* Fex, double* Phi_NObar, double* Phi_NHbar);
#endif
double Fex_H2O_FittedCorrelations::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const
{	double PhiEx = 0.0;
	//Quadratic part:
	ScalarFieldTilde V_O = double(gInfo.nr)*(COO*Ntilde[0] + COH*Ntilde[1]); Phi_Ntilde[0] += V_O;
	ScalarFieldTilde V_H = double(gInfo.nr)*(COH*Ntilde[0] + CHH*Ntilde[1]); Phi_Ntilde[1] += V_H;
	PhiEx += 0.5*gInfo.dV*(dot(V_O,Ntilde[0]) + dot(V_H,Ntilde[1]));

	//Compute gaussian weighted densities:
	ScalarField NObar = I(fex_gauss*Ntilde[0]), Phi_NObar; nullToZero(Phi_NObar, gInfo);
	ScalarField NHbar = I(fex_gauss*Ntilde[1]), Phi_NHbar; nullToZero(Phi_NHbar, gInfo);
	//Evaluated weighted density functional:
	#ifdef GPU_ENABLED
	ScalarField fex(ScalarFieldData::alloc(gInfo,isGpuEnabled()));
	Fex_H20_FittedCorrelations_gpu(gInfo.nr, NObar->dataGpu(), NHbar->dataGpu(),
		 fex->dataGpu(), Phi_NObar->dataGpu(), Phi_NHbar->dataGpu());
	PhiEx += integral(fex);
	#else
	PhiEx += gInfo.dV*threadedAccumulate(Fex_H2O_FittedCorrelations_calc, gInfo.nr,
		 NObar->data(), NHbar->data(), Phi_NObar->data(), Phi_NHbar->data());
	#endif
	//Convert gradients:
	Phi_Ntilde[0] += fex_gauss*Idag(Phi_NObar);
	Phi_Ntilde[1] += fex_gauss*Idag(Phi_NHbar);
	return PhiEx;
}

double Fex_H2O_FittedCorrelations::computeUniform(const double* N, double* Phi_N) const
{	return Fex_H2O_FittedCorrelations_calc(0, &N[0], &N[1], &Phi_N[0], &Phi_N[1]);
}
