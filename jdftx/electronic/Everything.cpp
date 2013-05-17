/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <electronic/Everything.h>
#include <electronic/ExactExchange.h>
#include <electronic/VanDerWaals.h>
#include <electronic/DOS.h>
#include <core/LatticeUtils.h>

void Everything::setup()
{
	//Removes unused pseudopotentials
	size_t counter = 0;
	while(counter < iInfo.species.size())
	{	if(iInfo.species[counter]->atpos.size() == 0)
		{	iInfo.species.erase(iInfo.species.begin()+counter);
		}
		else
			counter++;
	}
	
	//Symmetries (phase 1: lattice+basis dependent)
	symm.setup(*this);
	
	//Initialize the grid:
	gInfo.Gmax = sqrt(2*cntrl.Ecut); //Ecut = 0.5 Gmax^2
	gInfo.initialize(false, symm.getMatrices());

	//Exchange correlation setup
	logPrintf("\n---------- Exchange Correlation functional ----------\n");
	exCorr.setup(*this); //main functional
	if(exCorrDiff.size())
	{	logPrintf("\n---------- Auxiliary functionals for error estimation ----------\n");
		for(auto ec: exCorrDiff)
			ec->setup(*this); //comparison functionals evaluated at the end
	}

	//Atom positions, pseudopotentials
	iInfo.setup(*this);
	
	//Symmetries (phase 2: fftbox and k-mesh dependent)
	eInfo.kpointsFold();
	symm.setupMesh();
	
	//Set up k-points, bands and fillings
	eInfo.setup(*this, eVars.F, ener);

	//Set up the reduced bases for wavefunctions:
	logPrintf("\n----- Setting up reduced wavefunction bases (%s) -----\n",
		(cntrl.basisKdep==BasisKpointIndep) ? "single at Gamma point\n" :  "one per k-point");
	basis.resize(eInfo.nStates);
	double avg_nbasis = 0.0;
	for(int q=0; q<eInfo.nStates; q++)
	{	if(cntrl.basisKdep==BasisKpointDep)
			basis[q].setup(gInfo, iInfo, cntrl.Ecut, eInfo.qnums[q].k);
		else
		{	if(q==0) basis[q].setup(gInfo, iInfo, cntrl.Ecut, vector3<>(0,0,0));
			else basis[q] = basis[0];
		}
		avg_nbasis += 0.5*eInfo.qnums[q].weight * basis[q].nbasis;
	}
	logPrintf("average nbasis = %7.3lf , ideal nbasis = %7.3lf\n", avg_nbasis,
		pow(sqrt(2*cntrl.Ecut),3)*(gInfo.detR/(6*M_PI*M_PI)));
	logFlush();

	//Check if DOS calculator is needed:
	if(!dump.dos)
	{	for(auto dumpPair: dump)
			if(dumpPair.second == DumpDOS)
			{	dump.dos = std::make_shared<DOS>();
				break;
			}
	}
	
	//Check exact exchange requirements:
	if(exCorr.exxFactor()) //Check main functional first
		coulombParams.omegaSet.insert(exCorr.exxRange());
	for(auto ec: exCorrDiff) //Check the comparison functionals next
		if(ec->exxFactor())
			coulombParams.omegaSet.insert(ec->exxRange());
	if(dump.polarizability) coulombParams.omegaSet.insert(0.);
	
	//Coulomb-interaction setup (with knowledge of exact-exchange requirements):
	if(coulombParams.omegaSet.size() || dump.dos)
	{	//Initialize k-point sampled supercell:
		std::vector<vector3<>> kmeshUnreduced;
		for(const QuantumNumber& qnum: eInfo.qnums)
			kmeshUnreduced.push_back(qnum.k);
		coulombParams.supercell = std::make_shared<Supercell>(gInfo, kmeshUnreduced, symm.getMatrices(), symm.getKpointInvertList());
	}
	coulomb = coulombParams.createCoulomb(gInfo);
	
	//Exact exchange (if required)
	if(coulombParams.omegaSet.size())
		exx = std::make_shared<ExactExchange>(*this);

	//Setup VanDerWaals corrections
	if(iInfo.vdWenable || eVars.fluidParams.needsVDW())
		vanDerWaals = std::make_shared<VanDerWaals>(*this);
	
	//Setup wavefunctions, densities, fluid, output module etc:
	iInfo.update(ener); //needs to happen before eVars setup for LCAO
	eVars.setup(*this);
	dump.setup(*this);

	//Setup electronic minimization parameters:
	elecMinParams.nDim = 0;
	for(int s=0; s<eInfo.nStates; s++)
	{	elecMinParams.nDim += 2 * basis[s].nbasis * eInfo.nBands;
		if(eInfo.subspaceRotation)
			elecMinParams.nDim += eInfo.nBands * eInfo.nBands;
	}
	elecMinParams.fpLog = globalLog;
	elecMinParams.linePrefix = "ElecMinimize: ";
	elecMinParams.energyLabel = relevantFreeEnergyName(*this);
	
	//Setup ionic minimization parameters:
	ionicMinParams.nDim = 0;
	for(auto sp: iInfo.species)
		for(unsigned at=0; at<sp->atpos.size(); at++)
			ionicMinParams.nDim += sp->constraints[at].getDimension();
	if(!ionicMinParams.nDim) ionicMinParams.nDim = 1;
	ionicMinParams.fpLog = globalLog;
	ionicMinParams.linePrefix = "IonicMinimize: ";
	ionicMinParams.energyLabel = relevantFreeEnergyName(*this);
	
	//Setup fluid minimization parameters:
	switch(eVars.fluidParams.fluidType)
	{	case FluidLinearPCM:
		case FluidNonlocalPCM:
			fluidMinParams.nDim = gInfo.nr; break;
		case FluidNonlinearPCM: fluidMinParams.nDim = 4 * gInfo.nr; break;
		case FluidClassicalDFT:
			fluidMinParams.nDim = 4 * gInfo.nr; 
			break;
		default:
			fluidMinParams.nDim = 0;
	}
	fluidMinParams.fpLog = globalLog;
	fluidMinParams.linePrefix = "FluidMinimize: ";
	fluidMinParams.energyLabel = relevantFreeEnergyName(*this);
	//Tweak log format for fluid with inner minimization:
	if(eVars.fluidParams.fluidType==FluidLinearPCM
	|| eVars.fluidParams.fluidType==FluidNonlocalPCM)
	{	fluidMinParams.linePrefix = "\tFluidMinimize: ";
		if(!eVars.fluidParams.verboseLog)
			fluidMinParams.fpLog = nullLog;
	}
	
	//Setup lattice minimization parameters:
	latticeMinParams.fpLog = globalLog;
	latticeMinParams.linePrefix = "LatticeMinimize: ";
	latticeMinParams.energyLabel = relevantFreeEnergyName(*this);

	logPrintf("\n"); logFlush();
}
