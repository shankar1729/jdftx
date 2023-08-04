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

/**
 * @todo write docs
 */
class SpringConstant
{
public:
	bool calculateSpringConstant = false;
	std::vector<PerturbationInfo::Mode> modes;
	matrix kmatrix;
	
	SpringConstant(Everything& e);
	void setupModes();
	double perturbAtom(int atom);
	double computeMatrixElement(int a, int b);
	void getPerturbedEnergy(Energies& ener, PerturbationInfo::Mode modeA, PerturbationInfo::Mode modeB, double deltaA, double deltaB);
	double getIonDependentEnergy(const Energies& ener);
	void computeMatrix();
	
private:
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	IonInfo& iInfo;
	PerturbationInfo& pInfo;
	PerturbationSolver ps;
	
	int lastAtomPerturbed = -1;
	std::vector<ColumnBundle> Ctmp;
};

#endif // SPRINGCONSTANT_H
