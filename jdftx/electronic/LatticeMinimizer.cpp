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
	return *this;
}
void axpy(double alpha, const LatticeGradient& x, LatticeGradient& y)
{	y.lattice += alpha * x.lattice;
	axpy(alpha, x.ionic, y.ionic);
}
double dot(const matrix3<>& x, const matrix3<>& y)
{	return trace(x*(~y));
}
double dot(const LatticeGradient& x, const LatticeGradient& y)
{	return dot(x.lattice,y.lattice) + dot(x.ionic,y.ionic);
}
inline double nrm2(LatticeGradient& x) { return sqrt(dot(x,x)); }
LatticeGradient clone(const LatticeGradient& x) { return x; }
void randomize(LatticeGradient& x)
{	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			x.lattice(i,j) = Random::normal();
	randomize(x.ionic);
}

void bcast(matrix3<>& x)
{	for(int k=0; k<3; k++)
		mpiWorld->bcast(&x(k,0), 3);
}

//-------------  class LatticeMinimizer -----------------

LatticeMinimizer::LatticeMinimizer(Everything& e)
: e(e), imin(e), Rorig(e.gInfo.R), skipWfnsDrag(false)
{
	logPrintf("\n--------- Lattice Minimization ---------\n");
	
	//Ensure that lattice-move-scale is commensurate with symmetries:
	const std::vector<SpaceGroupOp>& sym = e.symm.getMatrices();
	for(const SpaceGroupOp& op: sym)
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(op.rot(i,j) && e.cntrl.lattMoveScale[i] != e.cntrl.lattMoveScale[j])
					die("latt-move-scale is not commensurate with symmetries:\n"
						"\t(Lattice vectors #%d and #%d are connected by symmetry,\n"
						"\tbut have different move scale factors %lg != %lg).\n",
						i, j, e.cntrl.lattMoveScale[i], e.cntrl.lattMoveScale[j]);
	
	//Check which lattice vectors can be altered:
	vector3<bool> isTruncated = e.coulombParams.isTruncated();
	for(int k=0; k<3; k++)
		isFixed[k] = (e.cntrl.lattMoveScale[k]==0.) || isTruncated[k];
}

void LatticeMinimizer::step(const LatticeGradient& dir, double alpha)
{	//Check if strain will become too large beforehand
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
	e.gInfo.R += alpha * dir.lattice;
	strain = e.gInfo.R * inv(Rorig) - matrix3<>(1.,1.,1.); //update strain to match current lattice vectors
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
	
	//! Compute energy (and ionic gradients if needed)
	imin.compute(grad ? &grad->ionic : 0, Kgrad ? &Kgrad->ionic : 0);
	
	//! Calculate lattice gradients (stresses) if necessary:
	if(grad)
	{	//Update IonInfo::stress (in Eh/a0^3 units):
		logPrintf("Calculating stress tensor ... "); logFlush();
		e.iInfo.computeStress();
		logPrintf(" done!\n");
		//Calculate grad->lattice (in Eh/a0 units):
		grad->lattice = (e.iInfo.stress * e.gInfo.detR) * e.gInfo.invRT;
		//Set Kgrad->lattice if necessary:
		if(Kgrad) Kgrad->lattice = grad->lattice;
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
{	logPrintf("\n"); e.gInfo.printLattice();
	logPrintf("\n# Strain Tensor:\n"); strain.print(globalLog, "%12lg ");
	logPrintf("\n# Stress Tensor:\n"); e.iInfo.stress.print(globalLog, "%12lg ");
	logPrintf("\n");
	return imin.report(iter); //IonicMinimizer::report will print atomic positions, forces etc.
}

void LatticeMinimizer::constrain(LatticeGradient& dir)
{	//Ionic part:
	imin.constrain(dir.ionic);
	/*
	//TODO: implement symmetrization and constraining of fixed directions
	
	*/
}

double LatticeMinimizer::safeStepSize(const LatticeGradient& dir) const
{	//Lattice criterion:
	matrix3<> dirStrain = dir.lattice * inv(Rorig);
	double alphaMaxLattice = 0.5 * GridInfo::maxAllowedStrain / nrm2(dirStrain);
	if(nrm2(strain) < GridInfo::maxAllowedStrain) //not already at a disallowed strain
		while(nrm2(strain+alphaMaxLattice*dirStrain) > GridInfo::maxAllowedStrain)
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
	e.coulomb = e.coulombParams.createCoulomb(e.gInfo, e.gInfoWfns, e.coulombWfns);
	e.iInfo.update(e.ener);
	logResume();
}
