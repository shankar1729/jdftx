/*-------------------------------------------------------------------
Copyright 2024 Ravishankar Sundararaman

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

#include <fluid/CANON.h>
#include <fluid/PCM_internal.h>
#include <core/ScalarFieldIO.h>
#include <core/SphericalHarmonics.h>
#include <electronic/Everything.h>


inline void setKernels(int i, double Gsq, double* preconditioner, double* metric)
{	preconditioner[i] = 1.0; 
	metric[i] = 1.0; //TODO: update for convergenece
}

CANON::CANON(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), Pulay(e.eVars.fluidParams.scfParams), NonlinearCommon(fsp, epsBulk)
{
	//Setup weight functions
	const double dG = gInfo.dGradial, Gmax = gInfo.GmaxGrid;
	unsigned nGradial = unsigned(ceil(Gmax/dG))+5;
	std::vector<double> w0_samples(nGradial), w1_samples(nGradial);
	for(unsigned i=0; i<nGradial; i++)
	{	double GR = (i*dG) * fsp.Res;
		double j0 = bessel_jl(0, GR);
		double j1_by_x_3 = j0 + bessel_jl(2, GR);  //3 j1(x)/x = j0(x) + j2(x)
		w0_samples[i] = fsp.Zcenter * (1.0 - j0);
        w1_samples[i] = j1_by_x_3;
	}
	w0.init(0, w0_samples, dG);
	w1.init(0, w1_samples, dG);

	//Setup the mixing and metric kernels
	preconditioner = std::make_shared<RealKernel>(gInfo);
	metric = std::make_shared<RealKernel>(gInfo);
	applyFuncGsq(gInfo, setKernels, preconditioner->data(), metric->data());
}

CANON::~CANON()
{	w0.free();
	w1.free();
}

void CANON::minimizeFluid()
{	logPrintf("\tCANON fluid occupying %lf of unit cell:", integral(shape[0])/gInfo.detR); logFlush();
	minimize();
}

void CANON::loadState(const char* filename)
{	nullToZero(phiTot, gInfo);
	loadRawBinary(phiTot, filename);
}

void CANON::saveState(const char* filename) const
{	saveRawBinary(phiTot, filename);
}

void CANON::dumpDensities(const char* filenamePattern) const
{	die("Not yet implemented.\n");
}

void CANON::set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	
	//Compute cavity shape function (0 to 1)
	nCavityNetTilde = nCavityTilde + getFullCore();
	nCavity = I(Sf[0] * nCavityNetTilde);
	updateCavity();

	//Initialize the state if it hasn't been loaded:
	if(!phiTot)
	{	phiTot = coulomb(rhoExplicitTilde);
		nullToZero(eps, gInfo);
		nullToZero(muPlus, gInfo);
		nullToZero(muMinus, gInfo);
	}
}


double CANON::get_Adiel_and_grad_internal(ScalarFieldTilde& grad_rhoExplicitTilde, ScalarFieldTilde& grad_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const
{	die("Not yet implemented.\n");
}


void CANON::getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const
{	die("Not yet implemented.\n");
}

double CANON::cycle(double dEprev, std::vector<double>& extraValues)
{
	//Get eps, mu from electrostatic potential:
	const VectorField Dphi = I(gradient(w1 * phiTot));
// 	callPref(dielectricEval->phiToState)(
// 		gInfo.nr, Dphi.dataPref(), shape[0]->dataPref(),
// 		gLookup, true, eps.dataPref(), NULL);
	double phiGzero = 0.0;
	if(screeningEval)
	{	//TODO: set phiGzero based on charge neutrality
		phiTot->setGzero(phiGzero);
		ScalarField phi = I(phiTot);
		callPref(screeningEval->phiToState)(
			gInfo.nr, phi->dataPref(), shape.back()->dataPref(),
			xLookup, true, muPlus->dataPref(), muMinus->dataPref(), NULL);
	}
	
	//Compute internal energy and bound charge:
	ScalarFieldTilde rhoFluidTilde; 
	if(screeningEval)
	{	//Compute ionic free energy and bound charge
		ScalarField Aout, rhoIon;
		initZero(Aout, gInfo);
		initZero(rhoIon, gInfo);
		callPref(screeningEval->freeEnergy)(gInfo.nr, 0.0, muPlus->dataPref(), muMinus->dataPref(), shape.back()->dataPref(),
			rhoIon->dataPref(), Aout->dataPref(), NULL, NULL, NULL);
		Adiel["Akappa"] = integral(Aout);
		rhoFluidTilde += J(rhoIon); //include bound charge due to ions
	}
	{	ScalarField Aout; initZero(Aout, gInfo);
		VectorField p; nullToZero(p, gInfo);
// 		callPref(dielectricEval->freeEnergy)(gInfo.nr, eps.const_dataPref(), shape[0]->dataPref(),
// 			p.dataPref(), Aout->dataPref(), vector3<double*>(NULL, NULL, NULL), NULL);
		Adiel["Aeps"] = integral(Aout);
		rhoFluidTilde -= w1 * divergence(J(p)); //include bound charge due to dielectric
	}
	
	//Compute new potential and electrostatic energies:
	ScalarFieldTilde phiFluidTilde = coulomb(rhoFluidTilde);
	ScalarFieldTilde phiExplicitTilde = coulomb(rhoExplicitTilde);
	Adiel["Coulomb"] = dot(rhoFluidTilde, O(0.5*phiFluidTilde + phiExplicitTilde));
	phiTot = phiExplicitTilde + phiFluidTilde;
	if(screeningEval) phiTot->setGzero(phiGzero);
	return Adiel;
}

void CANON::readVariable(ScalarFieldTilde& X, FILE* fp) const
{	nullToZero(X, gInfo);
	loadRawBinary(X, fp);
}

void CANON::writeVariable(const ScalarFieldTilde& X, FILE* fp) const
{	saveRawBinary(X, fp);
}

ScalarFieldTilde CANON::getVariable() const
{	return clone(phiTot);
}

void CANON::setVariable(const ScalarFieldTilde& X)
{	phiTot = clone(X);
}

ScalarFieldTilde CANON::precondition(const ScalarFieldTilde& X) const
{	return (*preconditioner) * X;
}

ScalarFieldTilde CANON::applyMetric(const ScalarFieldTilde& X) const
{	return (*metric) * X;
}
