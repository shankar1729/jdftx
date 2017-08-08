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
		mpiUtil->bcast(&x(k,0), 3);
}

//-------------  class LatticeMinimizer -----------------

LatticeMinimizer::LatticeMinimizer(Everything& e)
: e(e), imin(e), Rorig(e.gInfo.R), skipWfnsDrag(false)
{
	logPrintf("\n--------- Lattice Minimization ---------\n");
	
	//Ensure that lattice-move-scale is commensurate with symmetries:
	std::vector<SpaceGroupOp> sym = e.symm.getMatrices();
	for(const SpaceGroupOp& op: sym)
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(op.rot(i,j) && e.cntrl.lattMoveScale[i] != e.cntrl.lattMoveScale[j])
					die("latt-move-scale is not commensurate with symmetries:\n"
						"\t(Lattice vectors #%d and #%d are connected by symmetry,\n"
						"\tbut have different move scale factors %lg != %lg).\n",
						i, j, e.cntrl.lattMoveScale[i], e.cntrl.lattMoveScale[j]);
	
	//Check which lattice vectors can be altered:
	vector3<bool> isFixed, isTruncated = e.coulombParams.isTruncated();
	for(int k=0; k<3; k++)
		isFixed[k] = (e.cntrl.lattMoveScale[k]==0.) || isTruncated[k];
	
	//Create a orthonormal basis for strain commensurate with symmetries:
	for(int k=0; k<6; k++)
	{	//Initialize a basis element for arbitrary symmetric matrices:
		matrix3<int> s; //all zero:
		if(k<3) //diagonal strain
		{	s(k,k) = 1;
			if(isFixed[k]) continue; //strain alters fixed direction
		}
		else //off-diagonal strain
		{	int i=(k+1)%3;
			int j=(k+2)%3;
			s(i,j) = s(j,i) = 1;
			if(isFixed[i] || isFixed[j]) continue;  //strain alters fixed direction
		}
		//Symmetrize:
		matrix3<int> sSym;
		for(const SpaceGroupOp& op: sym)
		{	matrix3<int> mInv = det(op.rot) * adjugate(op.rot); //since |det(rot)| = 1
			sSym += mInv * s * op.rot;
		}
		//Orthonormalize w.r.t previous basis elements:
		matrix3<> strain(sSym); //convert from integer to double matrix
		for(const matrix3<>& sPrev: strainBasis)
			strain -= sPrev * dot(sPrev, strain);
		double strainNorm = nrm2(strain);
		if(strainNorm < symmThresholdSq) continue; //linearly dependent
		strainBasis.push_back((1./strainNorm) * strain);
	}
	if(!strainBasis.size())
		die("All lattice-vectors are constrained by coulomb truncation and/or\n"
			"latt-move-scale: please disable lattice minimization.\n");
	
	//Print initialization status:
	e.latticeMinParams.nDim = strainBasis.size();
	logPrintf("Minimization of dimension %lu over strains spanned by:\n", strainBasis.size());
	for(const matrix3<>& s: strainBasis)
	{	s.print(globalLog, " %lg ");
		logPrintf("\n");
	}

	h = 1e-5;
	
	//Set preconditioner (brings lattice Hessian to same dimensions as ionic):
	for(int iDir=0; iDir<3; iDir++)
		K[iDir] = e.cntrl.lattMoveScale[iDir] / e.gInfo.R.column(iDir).length();
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
	strain += alpha * dir.lattice;
	e.gInfo.R = Rorig + Rorig*strain; // Updates the lattice vectors to current strain
	bcast(e.gInfo.R); //ensure consistency to numerical precision
	updateLatticeDependent(e, true); // Updates lattice information but does not touch electronic state / calc electronic energy

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
		logPrintf("Calculating stress tensor... "); logFlush();
		calculateStress();
		logPrintf(" done!\n");
		//Calculate grad->lattice (in Eh units):
		grad->lattice = e.iInfo.stress * e.gInfo.detR;
		//Set Kgrad->lattice if necessary (in Eh/a0^2 units):
		if(Kgrad)
			Kgrad->lattice = Diag(K) * grad->lattice * Diag(K);
	}
	
	skipWfnsDrag = false; //computed at physical strain; safe to drag wfns at next step
	return relevantFreeEnergy(e);
}

void LatticeMinimizer::calculateStress()
{	matrix3<> E_strain;
	for(size_t i=0; i<strainBasis.size(); i++)
		E_strain += strainBasis[i]*centralDifference(strainBasis[i]);
	e.gInfo.R = Rorig + Rorig*strain;
	updateLatticeDependent(e);
	e.iInfo.stress = E_strain * (1./e.gInfo.detR);
}

double LatticeMinimizer::minimize(const MinimizeParams& params)
{	double result = Minimizable<LatticeGradient>::minimize(params);
	LatticeGradient dir; dir.ionic = e.iInfo.forces; //just needs to be right size; values irrelevant since used with step size 0
	step(dir, 0.); //so that population analysis may be performed at final positions / lattice
	return result;
}


double LatticeMinimizer::centralDifference(matrix3<> direction)
{ //! Implements a central difference derivative with O(h^4)
	
	e.gInfo.R = Rorig + Rorig*(strain+(-2*h*direction));
	updateLatticeDependent(e);
	const double En2h = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(strain+(-h*direction));
	updateLatticeDependent(e);
	const double Enh = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(strain+(h*direction));
	updateLatticeDependent(e);
	const double Eph = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(strain+(2*h*direction));
	updateLatticeDependent(e);
	const double Ep2h = relevantFreeEnergy(e);
	
	return (1./(12.*h))*(En2h - 8.*Enh + 8.*Eph - Ep2h);
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
	//Lattice part:
	matrix3<> latticeNew;
	for(const matrix3<>& s: strainBasis)
		latticeNew += s * dot(s, dir.lattice); //projection in basis
	dir.lattice = latticeNew;
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
{	mpiUtil->bcast(x);
	return x;
}

void LatticeMinimizer::updateLatticeDependent(Everything& e, bool ignoreElectronic)
{	logSuspend();
	e.gInfo.update();
	if(e.gInfoWfns)
	{	e.gInfoWfns->R = e.gInfo.R;
		e.gInfoWfns->update();
	}
	e.updateSupercell();
	e.coulomb = e.coulombParams.createCoulomb(e.gInfo);
	e.iInfo.update(e.ener);
	if(!ignoreElectronic)
	{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			e.eVars.orthonormalize(q);
		e.eVars.elecEnergyAndGrad(e.ener);
	}
	logResume();
}
