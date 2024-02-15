#include <fluid/CANON.h>
#include <core/ScalarFieldIO.h>
#include <electronic/Everything.h>


CANON::CANON(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), Pulay(e.eVars.fluidParams.scfParams)
{
}

CANON::~CANON()
{
}

void CANON::minimizeFluid()
{	minimize();
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

	//TODO: setup the mixing and metric kernels

	//Initialize the state if it hasn't been loaded:
	if(!phiTot) nullToZero(phiTot, gInfo);
}


double CANON::get_Adiel_and_grad_internal(ScalarFieldTilde& grad_rhoExplicitTilde, ScalarFieldTilde& grad_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const
{	die("Not yet implemented.\n");
}


void CANON::getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const
{	die("Not yet implemented.\n");
}

double CANON::cycle(double dEprev, std::vector<double>& extraValues)
{	die("Not yet implemented.\n");
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
