/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman

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

#include <wannier/DefectSupercell.h>
#include <wannier/Wannier.h>
#include <commands/parser.h>
#include <core/LatticeUtils.h>
#include <core/WignerSeitz.h>
#include <electronic/ColumnBundle.h>
#include <libgen.h>

void DefectSupercell::initialize(const Wannier* wannier)
{	this->wannier = wannier;
	this->e = wannier->e;
	logPrintf("\n\n########### Initializing supercell for defect '%s' ###########\n\n", name.c_str());

	//Change directory:
	//--- first store present directory
	const int bufSize = 1024;
	char currentDir[bufSize];
	getcwd(currentDir, bufSize);
	//--- split inFile into directory and filename within:
	char dirBuf[2][bufSize];
	for(int k=0; k<2; k++)
		strcpy(dirBuf[k], inFile.c_str());
	const char* inFileDir = dirname(dirBuf[0]);
	const char* inFileName = basename(dirBuf[1]);
	//--- switch to inFileDir:
	logPrintf("Switching to directory '%s'\n", inFileDir);
	if(chdir(inFileDir) != 0)
		die("Failed to switch to defect input directory '%s'\n\n", inFileDir);
	
	//Read input file:
	logPrintf("Reading supercell input file '%s' ... ", inFileName);
	std::vector< std::pair<string,string> > input = readInputFile(inFileName);
	logPrintf("done."); logFlush();
	eSup = std::make_shared<WannierEverything>(); //WannierEverything needed here to handle default wannier commands (though unused)
	parse(input, *eSup);
	
	//Alter relevant settings and initialize:
	//--- Remove unnecessary initializations:
	eSup->eVars.skipWfnsInit = true; //skip because wavefunctions are set from unit cell calculation
	eSup->eVars.wfnsFilename.clear();
	eSup->eVars.eigsFilename.clear();
	eSup->eVars.fluidInitialStateFilename.clear();
	eSup->eInfo.initialFillingsFilename.clear();
	eSup->scfParams.historyFilename.clear();
	eSup->eVars.initLCAO = false; //only unit cell state needed
	//--- Make sure electronic potentials are read in:
	eSup->cntrl.fixed_H = true;
	eSup->eVars.nFilenamePattern.clear();
	eSup->eVars.VFilenamePattern = eSup->dump.format;
	//--- Initialize (only need grid, pseudopotentials and Vscloc):
	eSup->Everything::setup();
	
	//Switch back to original working directory:
	logPrintf("Switching back to directory '%s'\n", currentDir);
	if(chdir(currentDir) != 0)
		die("Failed to switch back to unit cell input directory '%s'\n\n", currentDir);
	
	//Determine / check supercells:
	logPrintf("\n---- Defect supercell ----\n");
	//--- determine input supercell
	matrix3<> supInTmp = inv(e->gInfo.R) * eSup->gInfo.R;
	supIn = round(vector3<>(supInTmp(0,0), supInTmp(1,1), supInTmp(2,2)));
	logPrintf("Input supercell: "); supIn.print(globalLog, " %d ");
	double err = nrm2(supInTmp - Diag(vector3<>(supIn)));
	if(err > symmThreshold)
	{	logPrintf("Linear combinations of unit cell in defect supercell:\n");
		supInTmp.print(globalLog, " %lf ");
		die("Defect supercell is not diagonal/integral within threshold %lg\n\n", symmThreshold);
	}
	//--- check output supercell
	logPrintf("Output supercell: "); supOut.print(globalLog, " %d ");
	for(int j=0; j<3; j++)
		if(!supOut[j] or (e->eInfo.kfold[j] % supOut[j]))
		{	die("kpoint folding %d is not a multiple of output supercell count %d for lattice direction %d.\n",
				e->eInfo.kfold[j], supOut[j], j);
		}
	//--- map atom differences
	if(eSup->iInfo.coordsType == CoordsCartesian)
		xCenter = eSup->gInfo.invR * xCenter; //ensure lattice coordinates
	logPrintf("Center fractional coordinates: "); xCenter.print(globalLog, " %lf ");
	logPrintf("Atom changes from perfect to defective supercell:\n");
	logSuspend();
	WignerSeitz ws(eSup->gInfo.R);
	logResume();
	atposDefect.assign(eSup->iInfo.species.size(), std::vector<vector3<>>());
	atposRef = atposDefect;
	std::vector<bool> unitSpeciesFound(e->iInfo.species.size(), false);
	for(unsigned iSp=0; iSp<atposDefect.size(); iSp++)
	{	const SpeciesInfo& sp = *(eSup->iInfo.species[iSp]);
		//Map defect supercell atoms to Wigner-Seitz cell:
		for(vector3<> x: sp.atpos)
			atposDefect[iSp].push_back(xCenter + ws.restrict(x - xCenter));
		//Find corresponding unit cell species and atoms:
		unsigned nAdded = sp.atpos.size(); //number of atoms of this species added relative to perfect supercell
		unsigned nMoved = 0; //number of atoms of this species moved from perfect supercell
		double rSqMoved = 0.; //total square distance of atom movement
		double rCutMove = 1.; //threshold displacement in bohrs above which atom is considered added/deleted rather than moved
		for(unsigned iSpUnit=0; iSpUnit<e->iInfo.species.size(); iSpUnit++)
			if(sp.name == e->iInfo.species[iSpUnit]->name)
			{	unitSpeciesFound[iSpUnit] = true;
				//Add reference atom positions for perfect supercell:
				const SpeciesInfo& spUnit = *(e->iInfo.species[iSpUnit]);
				matrix3<> invSupIn = inv(Diag(vector3<>(supIn)));
				vector3<int> iR;
				for(iR[0]=0; iR[0]<supIn[0]; iR[0]++)
				for(iR[1]=0; iR[1]<supIn[1]; iR[1]++)
				for(iR[2]=0; iR[2]<supIn[2]; iR[2]++)
				{	for(vector3<> pos: spUnit.atpos)
					{	vector3<> x = invSupIn * (pos + iR);
						//Find corresponding closest atom in defect supercell:
						double rSqSmallest = DBL_MAX; unsigned iClosest = -1;
						for(unsigned iAtom=0; iAtom<sp.atpos.size(); iAtom++)
						{	vector3<> dx = x - atposDefect[iSp][iAtom];
							for(int dir=0; dir<3; dir++) dx[dir] -= floor(0.5 + dx[dir]); //minimum-image convention
							double rSq = eSup->gInfo.RTR.metric_length_squared(dx);
							if(rSq < rSqSmallest)
							{	rSqSmallest = rSq;
								iClosest = iAtom;
							}
						}
						//Count as a move if within threshold:
						if(rSqSmallest < std::pow(rCutMove,2))
						{	nAdded--;
							nMoved++;
							rSqMoved += rSqSmallest;
							const vector3<>& xRef = atposDefect[iSp][iClosest];
							x = xRef + ws.restrict(x - xRef);
						}
						else x = xCenter + ws.restrict(x - xCenter);
						atposRef[iSp].push_back(x);
					}
				}
				break;
			}
		unsigned nRemoved = atposRef[iSp].size() - nMoved;
		logPrintf("\t%s: added %u atoms, removed %u atoms, moved %u atoms",
			sp.name.c_str(), nAdded, nRemoved, nMoved);
		if(nMoved) logPrintf(" by %lf bohrs RMS", sqrt(rSqMoved/nMoved));
		logPrintf("\n");
		//Convert atoms to unit cell fratcional coordinates:
		matrix3<> DiagSupIn = Diag(vector3<>(supIn));
		for(vector3<>& x: atposDefect[iSp]) x = DiagSupIn * x;
		for(vector3<>& x: atposRef[iSp]) x = DiagSupIn * x;
	}
	for(unsigned iSpUnit=0; iSpUnit<e->iInfo.species.size(); iSpUnit++)
		if(not unitSpeciesFound[iSpUnit])
			logPrintf("WARNING: unit cell species %s not present in defect supercell.\n",
				e->iInfo.species[iSpUnit]->name.c_str());
		
	//TODO
}

