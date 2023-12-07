/*-------------------------------------------------------------------
Copyright 2022 Brandon Li

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

#ifndef PERTURB_PERTURBATIONSOLVER_H_
#define PERTURB_PERTURBATIONSOLVER_H_


#include <core/ScalarFieldArray.h>
#include <electronic/ColumnBundle.h>
#include <perturb/LinearSolver.h>

struct PerturbationGradient {
	std::vector<ColumnBundle> X; //!< First order change in unnormalized wfns
	const ElecInfo* eInfo;
	const class PerturbationInfo* pInfo;

	void init(const Everything& e);
	PerturbationGradient& operator*=(double alpha);
};

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y); //!< accumulate operation: Y += alpha*X
double dot(const PerturbationGradient& x, const PerturbationGradient& y); //!< inner product
double dot(const std::vector<ColumnBundle>& x, const std::vector<ColumnBundle>& y, PerturbationInfo* pInfo, ElecInfo *eInfo);
PerturbationGradient clone(const PerturbationGradient& x); //!< create a copy
void randomize(PerturbationGradient& x); //!< initialize with random numbers
ScalarField randomRealSpaceVector(ColumnBundle& ref, const  GridInfo* gInfo);

class PerturbationSolver : LinearSolvableIndefinite<PerturbationGradient>
{
public:
	PerturbationSolver(Everything& e);
	
	void hessian(PerturbationGradient& Av, const PerturbationGradient& v) override;
	void precondition(PerturbationGradient& v) override;
	void Ksqrt(PerturbationGradient& v) override;
	void solvePerturbation();
	
	void calcdGradTau(); //!< Get derivative of Grad E w.r.t. external perturbations (local potential, charge density, or atomic position)
	
	void getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y); //!< Get gradient of E w.r.t Y, used for finite difference tests.
	
	void computeIncommensurateWfns(); //!< Automatically set up and minimize incommensurate wfns, not implemented yet

	
	void getdn(ScalarFieldArray& dn, const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C = 0); //!< Derivative of density w.r.t. wavefunctions

	void getdnInc(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq); //!< Derivative of density w.r.t. wavefunctions, incommensurate case
	
	void getdnatom(ScalarFieldArray& dnatom); //!< Derivative of density w.r.t. atomic positions
	
	
    void updateExcorrCache(); //!< Update cache variables for calculation of LDA and GGA derivatives
	void updateHC(); //!< Update HC, OC, and density augmentation derivative. Run whever ionic positions change. 
	void updateNonlocalDerivs(); //!< Update Vnl derivatives, run when atom perturbation is applied.
	
	//! Apply Hamiltonian to arbitrary wavefunction, unlike in ElecVars::applyHamiltonian does not have to be eVars.C
	void applyH(const QuantumNumber& qnum, const diagMatrix& Fq, ColumnBundle& HCq, const ColumnBundle& Cq);
	
	//! Derivative of H operator w.r.t. Vscloc
	void dH_Vscloc(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq = 0);
	
	//! Derivative of H operator w.r.t. wavefunctions
	void dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq = 0);
	//! Derivative of H operator w.r.t. wavefunctions for incommensurate perturbations
	void dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const complexScalarFieldArray& dVscloc);
	
	//! Derivative of H operator w.r.t. external perturbations (local potential, charge density, or atomic position)
	void dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc);
	//! Derivative of H operator w.r.t. incommensurate external perturbations (local potential, charge density)
	void dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, int qsign);
	
	//! Derivative of Vscloc w.r.t density
	void getdVsclocPsi(const ScalarFieldArray dn, ScalarFieldArray& dVscloc);
	//! Derivative of Vscloc w.r.t density, incommensurate
	void getdVsclocPsi(const complexScalarFieldArray dn, complexScalarFieldArray& dVscloc, vector3<>* q);
	//! Derivative of Vscloc w.r.t external perturbations
	void getdVsclocTau(ScalarFieldArray& dVscloc, ScalarFieldArray* dn = 0);
	//! Derivative of Vscloc w.r.t external perturbations, incommensurate
	void getdVsclocTau(complexScalarFieldArray& dVscloc, int qsign);
	
	//! Derivative of nXC w.r.t core density
	ScalarFieldArray getdnXC(const ScalarField dnCore) const;
	
	//! Copy of ElecVars::get_nTot for any n
	ScalarField get_nTot(const ScalarFieldArray n) const { return n.size()==1 ? clone(n[0]) : n[0]+n[1]; }
	complexScalarField get_nTot(const complexScalarFieldArray n) const { return n.size()==1 ? clone(n[0]) : n[0]+n[1]; }

private:
	Everything& e;
	class ElecVars& eVars;
	ElecInfo& eInfo;
	IonInfo& iInfo;
	PerturbationInfo& pInfo;
};

void init(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const Basis* basis, const ElecInfo* eInfo);

#endif /* PERTURB_PERTURBATIONSOLVER_H_ */
