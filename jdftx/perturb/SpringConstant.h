/*-------------------------------------------------------------------
Copyright 2023 Brandon Li

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

#ifndef PERTURB_SPRINGCONSTANT_H
#define PERTURB_SPRINGCONSTANT_H

#include <perturb/PerturbationSolver.h>
#include <perturb/PerturbationInfo.h>
#include <electronic/IonicMinimizer.h>

class SpringConstant
{
public:
	SpringConstant(Everything& e);
	double computeMatrixElement(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Compute individual matrix element
	
	double dsqEpair(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of pair potentials (Ewald and VDW)
	double dsqEnl(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of nonlocal energy
	double dsqEloc(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of local energy
	
	void getPerturbedEnergy(Energies& ener, std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB, double deltaA, double deltaB); //!< Perturb atoms by amounts deltaA and deltaB in the directions determined by modeA and modeB respectively, then compute energy
	IonicGradient getPhononMatrixColumn(std::shared_ptr<AtomPerturbation> mode, double dr = 0.0); //!< Compute part of OmegaSq corresponding to single atom perturbation, used by phonon package
	
private:
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	IonInfo& iInfo;
	PerturbationInfo& pInfo;
	PerturbationSolver ps;
	
	std::vector<ColumnBundle> Ctmp; //!< Temporary storage for eVars.C during phonon calculation
};

#endif // PERTURB_SPRINGCONSTANT_H
