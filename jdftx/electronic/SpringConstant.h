/*
 * Copyright 2023 <copyright holder> <email>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SPRINGCONSTANT_H
#define SPRINGCONSTANT_H

#include <electronic/PerturbationSolver.h>
#include <electronic/PerturbationInfo.h>
#include <electronic/IonicMinimizer.h>

/**
 * @todo write docs
 */
class SpringConstant
{
public:
	bool calculateSpringConstant = false;
	std::vector<std::shared_ptr<AtomPerturbation>> modes;
	matrix kmatrix;
	
	SpringConstant(Everything& e);
	void setupModes();
	double perturbAtom(int atom);
	double computeMatrixElement(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	
	
	std::vector<ColumnBundle> dsqC;
	std::vector<matrix> dsqVdagC;
	ScalarFieldArray dsqn;
	ScalarField dsqE_naug;
	double dsqQuantities(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	
	double dsqEpair(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	double dsqEnl(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	double dsqEloc(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	double dsqExc(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	double dsqExccore(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	double dsqEH(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB);
	
	void getPerturbedEnergy(Energies& ener, std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB, double deltaA, double deltaB);
	void computeSubMatrix();
	IonicGradient getPhononMatrixColumn(std::shared_ptr<AtomPerturbation> mode, double dr = 0.0);
	
private:
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	IonInfo& iInfo;
	PerturbationInfo& pInfo;
	PerturbationSolver ps;
	
	std::vector<ColumnBundle> Ctmp;
};

#endif // SPRINGCONSTANT_H
