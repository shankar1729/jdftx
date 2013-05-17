/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#include <fluid/Fex_ScalarEOS_internal.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_LJ.h>
#include <core/Units.h>
#include <core/Operators.h>
#include <electronic/operators.h>

string rigidMoleculeCDFT_ScalarEOSpaper = "R. Sundararaman and T.A. Arias, arXiv:1302.0026";

Fex_ScalarEOS::Fex_ScalarEOS(const FluidMixture* fluidMixture, const FluidComponent* comp, const ScalarEOS& eos)
: Fex(fluidMixture, comp), eos(eos)
{
	double Rhs = 0.;
	for(const auto& s: molecule.sites)
		if(s->Rhs)
		{	assert(!Rhs); //Ensure single hard sphere site
			Rhs = s->Rhs;
		}
	assert(Rhs);
	Vhs = (4*M_PI/3) * pow(Rhs,3);

	//Initialize the mean field kernel:
	double sigma = 2.*Rhs;
	setLJatt(fex_LJatt, gInfo, -9.0/(32*sqrt(2)*M_PI*pow(sigma,3)), sigma);
	Citations::add("Scalar-EOS liquid functional", rigidMoleculeCDFT_ScalarEOSpaper);
}
Fex_ScalarEOS::~Fex_ScalarEOS()
{	fex_LJatt.free();
}

double Fex_ScalarEOS::compute(const DataGptr* Ntilde, DataGptr* Phi_Ntilde) const
{	//Compute LJatt weighted density:
	DataRptr Nbar = I(fex_LJatt*Ntilde[0]);
	//Evaluated weighted density functional:
	DataRptr Aex, Aex_Nbar; nullToZero(Aex, gInfo); nullToZero(Aex_Nbar, gInfo);
	callPref(eos.evaluate)(gInfo.nr, Nbar->dataPref(), Aex->dataPref(), Aex_Nbar->dataPref(), Vhs);
	//Convert gradients:
	DataRptr NO = I(Ntilde[0]);
	Phi_Ntilde[0] += fex_LJatt*Idag(NO*Aex_Nbar) + Idag(Aex);
	return gInfo.dV*dot(NO,Aex);
}

double Fex_ScalarEOS::computeUniform(const double* N, double* Phi_N) const
{	double AexPrime, Aex;
	eos.evaluate(1, &N[0], &Aex, &AexPrime, Vhs);
	Phi_N[0] += Aex + N[0]*AexPrime;
	return N[0]*Aex;
}


//--------------- class JeffereyAustinEOS ---------------

JeffereyAustinEOS::JeffereyAustinEOS(double T)
: eval(std::make_shared<JeffereyAustinEOS_eval>(T))
{
}

double JeffereyAustinEOS::vdwRadius() const
{	return eval->vdwRadius();
}

void evalJeffereyAustinEOS_sub(size_t iStart, size_t iStop, const double* N, double* Aex, double* Aex_N, double Vhs, const JeffereyAustinEOS_eval& eval)
{	for(size_t i=iStart; i<iStop; i++) eval(i, N, Aex, Aex_N, Vhs);
}
void JeffereyAustinEOS::evaluate(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const
{	threadLaunch(evalJeffereyAustinEOS_sub, nData, N, Aex, Aex_N, Vhs, *eval);
}
#ifdef GPU_ENABLED
void evalJeffereyAustinEOS_gpu(int nr, const double* Nbar, double* Fex, double* Phi_Nbar, double Vhs, const JeffereyAustinEOS_eval& eval);
void JeffereyAustinEOS::evaluate_gpu(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const
{	evalJeffereyAustinEOS_gpu(nData, N, Aex, Aex_N, Vhs, *eval);
}
#endif

//--------------- class TaoMasonEOS ---------------

TaoMasonEOS::TaoMasonEOS(double T, double Tc, double Pc, double omega)
: eval(std::make_shared<TaoMasonEOS_eval>(T, Tc, Pc, omega))
{
}

double TaoMasonEOS::vdwRadius() const
{	return eval->vdwRadius();
}

void evalTaoMasonEOS_sub(size_t iStart, size_t iStop, const double* N, double* Aex, double* Aex_N, double Vhs, const TaoMasonEOS_eval& eval)
{	for(size_t i=iStart; i<iStop; i++) eval(i, N, Aex, Aex_N, Vhs);
}
void TaoMasonEOS::evaluate(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const
{	threadLaunch(evalTaoMasonEOS_sub, nData, N, Aex, Aex_N, Vhs, *eval);
}
#ifdef GPU_ENABLED
void evalTaoMasonEOS_gpu(int nr, const double* Nbar, double* Fex, double* Phi_Nbar, double Vhs, const TaoMasonEOS_eval& eval);
void TaoMasonEOS::evaluate_gpu(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const
{	evalTaoMasonEOS_gpu(nData, N, Aex, Aex_N, Vhs, *eval);
}
#endif
