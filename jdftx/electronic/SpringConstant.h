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

#ifndef SPRINGCONSTANT_H
#define SPRINGCONSTANT_H

#include <electronic/PerturbationSolver.h>
#include <electronic/PerturbationInfo.h>
#include <electronic/IonicMinimizer.h>

class SpringConstant
{
public:
	bool calculateSpringConstant = false; //!< True if at least one of the atoms specified by the ion command has the 'spring' flag enabled
	std::vector<std::shared_ptr<AtomPerturbation>> modes; //!< List of perturbations with three DOF for each perturbed atom
	matrix kmatrix; //!< Matrix with entries containing second derivative of total energy w.r.t. pairs of modes
	
	SpringConstant(Everything& e);
	void setupModes(); //!< Run after ion commands
	double computeMatrixElement(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Compute individual matrix element

	//double dsqQuantities(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	
	double dsqEpair(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of pair potentials (Ewald and VDW)
	double dsqEnl(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of nonlocal energy
	double dsqEloc(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB); //!< Second derivative of local energy
	
	void getPerturbedEnergy(Energies& ener, std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB, double deltaA, double deltaB); //!< Perturb atoms by amounts deltaA and deltaB in the directions determined by modeA and modeB respectively, then compute energy
	void computeSubMatrix(); //!< Compute a submatrix of the spring constant matrix determined by SpringConstat::modes
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

#endif // SPRINGCONSTANT_H
