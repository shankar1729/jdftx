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

#include <fluid/Fex_H2O_FittedCorrelations_internal.h>
#include <fluid1D/Fex_H2O_FittedCorrelations.h>
#include <fluid1D/FluidMixture.h>
#include <core1D/Operators.h>
#include <core/Units.h>
#include <core/Util.h>

static const double rOH = 1.0*Angstrom;
static const double thetaHOH = acos(-1.0/3);

inline void setKernels(int i, const double* Garr, double* COO, double* COH, double* CHH, double* fex_gauss, double* siteChargeKernel)
{
	//Correlation kernels (gaussians of magnitude *A[:] centered at *B[:] with width *C[:]):
	static const double COO_A[6] = {-0.0271153, -0.0795576, 0.096648, -0.0291517, 0.0227052, -0.0109078};
	static const double COO_B[6] = {1.0569, 1.63853, 1.87815, 2.46787, 3.10047, 3.7162};
	static const double COO_C[6] = {0.0204818, 0.129513, 0.0736776, 0.0563636, 0.0844132, 0.0452023};

	static const double COH_A[6] = {-0.00801259, 0.0379836, -0.0380737, 0.0254576, -0.00291329, 0.00109967};
	static const double COH_B[6] = {1.19859, 1.72897, 2.15839, 2.72112, 3.29894, 3.7659};
	static const double COH_C[6] = {0.028028, 0.117534, 0.213336, 0.154157, 0.0265282, 0.0315518};

	static const double CHH_A[2] = {-0.013959, 0.0295776};
	static const double CHH_B[2] = {1.88697, 2.53164};
	static const double CHH_C[2] = {0.104186, 0.0869848};

	const double Gc = 0.33; //parameter for crossover of coulomb matrix
	const double r0 = 4.2027; //width of excess function correlations

	double G = Garr[i];

	//The electrostatic rolloff from the paper becoems a site charge profile here:
	siteChargeKernel[i] = 1.0/sqrt(1 + pow(G/Gc,4));

	//Fitted correlations:
	COO[i] = 0.0; for(int j=0; j<6; j++) COO[i] += COO_A[j] * exp(-pow(G-COO_B[j], 2) / COO_C[j]);
	COH[i] = 0.0; for(int j=0; j<6; j++) COH[i] += COH_A[j] * exp(-pow(G-COH_B[j], 2) / COH_C[j]);
	CHH[i] = 0.0; for(int j=0; j<2; j++) CHH[i] += CHH_A[j] * exp(-pow(G-CHH_B[j], 2) / CHH_C[j]);

	//convolution kernel for the weighted density:
	fex_gauss[i] = exp(-pow(0.5*r0*G, 2));
}


Fex_H2O_FittedCorrelations::Fex_H2O_FittedCorrelations(FluidMixture& fluidMixture)
: Fex(fluidMixture),
COO(gInfo.S), COH(gInfo.S), CHH(gInfo.S), fex_gauss(gInfo.S), siteChargeKernel(gInfo.S),
propO(gInfo, 0.0,0.0, 0.8476,&siteChargeKernel),
propH(gInfo, 0.0,0.0,-0.4238,&siteChargeKernel),
molecule("H2O",
	&propO,
		 vector3<>(0,0,0),
	&propH,
		 vector3<>(0, -rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)),
		 vector3<>(0, +rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)) )
{
	if(fabs(T/Kelvin-298)>1)
		die("The FittedCorrelations functional is only valid at T=298K.\n")
	
	//Initialize the kernels:
	serialLoop(setKernels, gInfo.S, gInfo.G.data(),
		COO.data(), COH.data(), CHH.data(), fex_gauss.data(), siteChargeKernel.data());
}

double Fex_H2O_FittedCorrelations::get_aDiel() const
{	return 0.959572098592; //Evaluating eps/(eps-1) - epsNI/(epsNI-1) for eps=78.4 and epsNI at STP
}

//1D wrapper to the function in the 3D code - multiplies by the non-uniform quadrature weights
inline double Fex_H2O_FittedCorrelations_calc1D(int i, const double* wArr, const double* NObar, const double* NHbar, double* grad_NObar, double* grad_NHbar)
{	double Fex = Fex_H2O_FittedCorrelations_calc(i, NObar, NHbar, grad_NObar, grad_NHbar);
	double w = wArr[i];
	grad_NObar[i] *= w;
	grad_NHbar[i] *= w;
	return Fex * w;
}

double Fex_H2O_FittedCorrelations::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const
{	double PhiEx = 0.0;
	//Quadratic part:
	ScalarFieldTilde V_O = O(COO*Ntilde[0] + COH*Ntilde[1]); grad_Ntilde[0] += V_O;
	ScalarFieldTilde V_H = O(COH*Ntilde[0] + CHH*Ntilde[1]); grad_Ntilde[1] += V_H;
	PhiEx += 0.5*(dot(V_O,Ntilde[0]) + dot(V_H,Ntilde[1]));

	//Compute gaussian weighted densities:
	ScalarField NObar = I(fex_gauss*Ntilde[0]), grad_NObar; nullToZero(grad_NObar, gInfo);
	ScalarField NHbar = I(fex_gauss*Ntilde[1]), grad_NHbar; nullToZero(grad_NHbar, gInfo);
	//Evaluated weighted density functional:
	PhiEx += serialAccumulate(Fex_H2O_FittedCorrelations_calc1D, gInfo.S, gInfo.w.data(),
		 NObar.data(), NHbar.data(), grad_NObar.data(), grad_NHbar.data());
	//Convert gradients:
	grad_Ntilde[0] += fex_gauss*Idag(grad_NObar);
	grad_Ntilde[1] += fex_gauss*Idag(grad_NHbar);
	return PhiEx;
}

double Fex_H2O_FittedCorrelations::computeUniform(const double* N, double* grad_N) const
{	return Fex_H2O_FittedCorrelations_calc(0, &N[0], &N[1], &grad_N[0], &grad_N[1]);
}

void Fex_H2O_FittedCorrelations::directCorrelations(const double* N, ScalarFieldTildeCollection& C) const
{	//Compute second derivatives of free energy density:
	double NO = N[0], NH = N[1], dN = 1e-7*NO;
	NO += 1*dN; double Op_aO, Op_aH; Fex_H2O_FittedCorrelations_calc(0, &NO, &NH, &Op_aO, &Op_aH);
	NO -= 2*dN; double Om_aO, Om_aH; Fex_H2O_FittedCorrelations_calc(0, &NO, &NH, &Om_aO, &Om_aH);
	NO += 1*dN;
	NH += 1*dN; double Hp_aO, Hp_aH; Fex_H2O_FittedCorrelations_calc(0, &NO, &NH, &Hp_aO, &Hp_aH);
	NH -= 2*dN; double Hm_aO, Hm_aH; Fex_H2O_FittedCorrelations_calc(0, &NO, &NH, &Hm_aO, &Hm_aH);
	NH += 1*dN;
	double a_OO = (Op_aO - Om_aO) / (2*dN);
	double a_OH = (Op_aH - Om_aH + Hp_aO - Hm_aO) / (4*dN);
	double a_HH = (Hp_aH - Hm_aH) / (2*dN);
	//Accumulate correlations:
	ScalarFieldTilde fex_gaussSq = fex_gauss * ScalarFieldTilde(fex_gauss, gInfo); //square of weight function
	C[fluidMixture.corrFuncIndex(0,0,this)] += (ScalarFieldTilde(COO,gInfo) + a_OO * fex_gaussSq);
	C[fluidMixture.corrFuncIndex(0,1,this)] += (ScalarFieldTilde(COH,gInfo) + a_OH * fex_gaussSq);
	C[fluidMixture.corrFuncIndex(1,1,this)] += (ScalarFieldTilde(CHH,gInfo) + a_HH * fex_gaussSq);
}
