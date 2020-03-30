/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Ravishankar Sundararaman

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

#include <electronic/LatticeMinimizer.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/LatticeUtils.h>
#include <core/Random.h>

//Functions required by Minimizable<LatticeGradient>
LatticeGradient& LatticeGradient::operator*=(double scale)
{	lattice *= scale;
	ionic *= scale;
	for(double& v: thermo) v *= scale;
	return *this;
}

LatticeGradient& LatticeGradient::operator+=(const LatticeGradient& other)
{	axpy(1., other, *this);
	return *this;
}

LatticeGradient LatticeGradient::operator+(const LatticeGradient& other) const
{	LatticeGradient result(*this);
	axpy(1., other, result);
	return result;
}

LatticeGradient LatticeGradient::operator-(const LatticeGradient& other) const
{	LatticeGradient result(*this);
	axpy(-1., other, result);
	return result;
}

void axpy(double alpha, const LatticeGradient& x, LatticeGradient& y)
{	y.lattice += alpha * x.lattice;
	axpy(alpha, x.ionic, y.ionic);
	//Optional thermostat DOFs used only by IonicDynamics:
	assert(x.thermo.size() == y.thermo.size());
	for(size_t i=0; i<x.thermo.size(); i++)
		y.thermo[i] += alpha * x.thermo[i];
}
double dot(const matrix3<>& x, const matrix3<>& y)
{	return trace(x*(~y));
}
double dot(const LatticeGradient& x, const LatticeGradient& y)
{	double result = dot(x.lattice,y.lattice) + dot(x.ionic,y.ionic);
	//Optional thermostat DOFs used only by IonicDynamics:
	assert(x.thermo.size() == y.thermo.size());
	for(size_t i=0; i<x.thermo.size(); i++)
		result += x.thermo[i] * y.thermo[i];
	return result;
}
inline double nrm2(LatticeGradient& x) { return sqrt(dot(x,x)); }
LatticeGradient clone(const LatticeGradient& x) { return x; }
void randomize(LatticeGradient& x)
{	for(int i=0; i<3; i++)
		for(int j=i; j<3; j++)
			x.lattice(i,j) = (x.lattice(j,i) = Random::normal());
	randomize(x.ionic);
	//Optional thermostat DOFs used only by IonicDynamics:
	for(double& v: x.thermo) v = Random::normal();
}

void bcast(matrix3<>& x)
{	for(int k=0; k<3; k++)
		mpiWorld->bcast(&x(k,0), 3);
}

//-------------  class LatticeMinimizer -----------------

LatticeMinimizer::LatticeMinimizer(Everything& e, bool dynamicsMode, bool statP, bool statStress)
: e(e), dynamicsMode(dynamicsMode), statP(statP), statStress(statStress),
	imin(e, dynamicsMode), Rorig(e.gInfo.R), skipWfnsDrag(false)
{
	if(not dynamicsMode) logPrintf("\n--------- Lattice Minimization ---------\n");
	
	//Ensure that lattice-move-scale is commensurate with symmetries:
	const std::vector<SpaceGroupOp>& sym = e.symm.getMatrices();
	for(const SpaceGroupOp& op: sym)
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(op.rot(i,j) && e.cntrl.lattMoveScale[i] != e.cntrl.lattMoveScale[j])
					die("latt-move-scale is not commensurate with symmetries:\n"
						"\tLattice vectors #%d and #%d are connected by symmetry\n"
						"\tbut have different move scale factors %lg != %lg.\n\n",
						i+1, j+1, e.cntrl.lattMoveScale[i], e.cntrl.lattMoveScale[j]);
	
	//Check which lattice vectors can be altered:
	const vector3<bool> isTruncated = e.coulombParams.isTruncated();
	std::vector<vector3<>> vFixed; //unit vectors spanning fixed subspace
	Pfree = matrix3<>(1,1,1); //initially identity (no projection)
	double lattMoveScaleMean = 0.; //mean value of free directions
	for(int k=0; k<3; k++)
		if((e.cntrl.lattMoveScale[k]==0.) || isTruncated[k])
		{	vector3<> vk = Rorig.column(k);
			for(vector3<> vPrev: vFixed)
				vk -= vPrev*dot(vPrev, vk); //orthogonal to previous fixed vectors
			vk *= (1./vk.length()); //normalize
			vFixed.push_back(vk); //remember this direction to orthogonalize future ones
			Pfree -= outer(vk, vk); //project out this direction
		}
		else lattMoveScaleMean += e.cntrl.lattMoveScale[k];
	if(vFixed.size()==3 and (not dynamicsMode))
		die("No lattice directions free for lattice minimization due to truncation and/or latt-move-scale.\n\n");
	lattMoveScaleMean /= (3.-vFixed.size());
	
	//Preconditioning factor:
	latticeK = std::pow(lattMoveScaleMean / std::pow(fabs(det(Rorig)),1./3), 2); //match dimensions to ionic Hessian
}

void LatticeMinimizer::step(const LatticeGradient& dir, double alpha)
{	if(dynamicsMode and (not (statP or statStress)))
	{	imin.step(dir.ionic, alpha); //since lattice constant, bypass more expensive processing below
		return;
	}
	
	//Check if strain will become too large beforehand
	//(so that we can avoid updating wavefunctions for steps that will fail)
	if(nrm2(strain+alpha*dir.lattice) > GridInfo::maxAllowedStrain) //strain will become large
		skipWfnsDrag = true; //skip wavefunction drag till a 'real' compute occurs at an acceptable strain
	
	//Update atomic positions first (with associated wavefunction drag, if any):
	imin.step(dir.ionic, alpha);
	
	//Project wavefunctions to atomic orbitals:
	std::vector<matrix> coeff(e.eInfo.nStates); //best fit coefficients
	int nAtomic = e.iInfo.nAtomicOrbitals();
	if(e.cntrl.dragWavefunctions && nAtomic && (!skipWfnsDrag))
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	//Get atomic orbitals for old lattice:
			ColumnBundle psi = e.iInfo.getAtomicOrbitals(q, false);
			//Fit the wavefunctions to atomic orbitals (minimize C0^OC0 where C0 is the remainder)
			ColumnBundle Opsi = O(psi); //non-trivial cost for uspp
			coeff[q] = inv(psi^Opsi) * (Opsi^e.eVars.C[q]);
			Opsi.free();
			e.eVars.C[q] -= psi * coeff[q]; //now contains residual
		}
	
	//Change lattice:
	const matrix3<> id(1,1,1); //identity
	matrix3<> strainFactor = id + alpha*dir.lattice;
	e.gInfo.R = strainFactor * e.gInfo.R;
	strain = e.gInfo.R * inv(Rorig) - id;
	bcast(e.gInfo.R); //ensure consistency to numerical precision
	bcast(strain); //ensure consistency to numerical precision
	updateLatticeDependent(e); // Updates lattice information

	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Restore wavefunctions from atomic orbitals:
		if(e.cntrl.dragWavefunctions && nAtomic && (!skipWfnsDrag))
		{	//Get atomic orbitals for new lattice:
			ColumnBundle psi = e.iInfo.getAtomicOrbitals(q, false);
			//Reconstitute wavefunctions:
			e.eVars.C[q] += psi * coeff[q];
		}
		e.eVars.orthonormalize(q); //Reorthonormalize wavefunctions
	}
}

double LatticeMinimizer::compute(LatticeGradient* grad, LatticeGradient* Kgrad)
{
	//Check for large lattice strain
	if(nrm2(strain) > GridInfo::maxAllowedStrain)
	{	logPrintf("\nBacking off lattice step because strain tensor has become enormous:\n"); strain.print(globalLog, "%10lg ");
		logPrintf("If such large strain is expected, restart calculation with these lattice vectors to prevent Pulay errors:\n");
		e.gInfo.printLattice();
		logPrintf("\n");
		return NAN;
	}
	if(not e.iInfo.checkPositions())
	{	logPrintf("\nBacking off lattice step since it caused pseudopotential core overlaps.\n");
		return NAN;
	}
	
	//! Compute energy (and ionic gradients, stress if needed)
	imin.compute(grad ? &grad->ionic : 0, Kgrad ? &Kgrad->ionic : 0);
	
	//! Calculate lattice gradients (from stress computed along with forces above) if necessary:
	if(grad)
	{	//Calculate grad->lattice (in Eh units):
		grad->lattice = e.iInfo.stress * e.gInfo.detR;
		//Set Kgrad->lattice if necessary:
		if(Kgrad)
		{	Kgrad->lattice = grad->lattice * latticeK;
			constrain(*Kgrad);
		}
	}
	
	skipWfnsDrag = false; //computed at physical strain; safe to drag wfns at next step
	return relevantFreeEnergy(e);
}

double LatticeMinimizer::minimize(const MinimizeParams& params)
{	double result = Minimizable<LatticeGradient>::minimize(params);
	LatticeGradient dir; dir.ionic = e.iInfo.forces; //just needs to be right size; values irrelevant since used with step size 0
	step(dir, 0.); //so that population analysis may be performed at final positions / lattice
	return result;
}

bool LatticeMinimizer::report(int iter)
{	if((not dynamicsMode) or statP or statStress)
	{	logPrintf("# Lattice vectors:\n"); e.gInfo.printLattice();
		logPrintf("\n# Strain tensor in Cartesian coordinates:\n"); strain.print(globalLog, "%12lg ", true, 1e-14);
	}
	return imin.report(iter); //IonicMinimizer::report will print stress, atomic positions, forces etc.
}

void LatticeMinimizer::constrain(LatticeGradient& dir)
{	//Ionic part:
	imin.constrain(dir.ionic);
	//Lattice part:
	dir.lattice = 0.5*(dir.lattice + ~(dir.lattice)); //ensure symmetric tensor
	dir.lattice = Pfree * dir.lattice * Pfree; //latt-move-scale and truncation constraints
	e.symm.symmetrize(dir.lattice); //lattice symmetries
}

double LatticeMinimizer::safeStepSize(const LatticeGradient& dir) const
{	//Lattice criterion:
	double alphaMaxLattice = 0.5 * GridInfo::maxAllowedStrain / nrm2(dir.lattice);
	if(nrm2(strain) < GridInfo::maxAllowedStrain) //not already at a disallowed strain
		while(nrm2(strain+alphaMaxLattice*dir.lattice) > GridInfo::maxAllowedStrain)
			alphaMaxLattice *= 0.5; //reduce step size further till new position will become safe
	//Ionic criterion:
	double alphaMaxIonic = imin.safeStepSize(dir.ionic);
	return std::min(alphaMaxLattice, alphaMaxIonic); //minimum of the two criterio
}

double LatticeMinimizer::sync(double x) const
{	mpiWorld->bcast(x);
	return x;
}

void LatticeMinimizer::updateLatticeDependent(Everything& e)
{	logSuspend();
	e.gInfo.update();
	if(e.gInfoWfns)
	{	e.gInfoWfns->R = e.gInfo.R;
		e.gInfoWfns->update();
	}
	e.updateSupercell();
	e.coulombParams.recreateCoulomb(e.gInfo, e.gInfoWfns, e.coulomb, e.coulombWfns);
	e.iInfo.update(e.ener);
	logResume();
}
