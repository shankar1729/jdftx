/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <wannier/Wannier.h>
#include <wannier/WannierMinimizerFD.h>
#include <wannier/WannierMinimizerRS.h>

Wannier::Wannier() : needAtomicOrbitals(false), localizationMeasure(LM_FiniteDifference),
	bStart(0), outerWindow(false), innerWindow(false), mainWindow(false), nMain(0),
	saveWfns(false), saveWfnsRealSpace(false), saveMomenta(false),
	loadRotations(false), numericalOrbitalsOffset(0.5,0.5,0.5)
{
}

void Wannier::setup(const Everything& everything)
{	e = &everything;
	logPrintf("\n---------- Initializing Wannier Function solver ----------\n");
	//Initialize minimization parameters:
	minParams.fpLog = globalLog;
	minParams.linePrefix = "WannierMinimize: ";
	minParams.energyLabel = "Omega";
	//Check window settings:
	if(mainWindow)
	{	if(!innerWindow) die("Main window requires that an inner window be specified.\n");
		if(eMainMin<eInnerMin || eMainMax>eInnerMax)
			die("Main window must lie entirely within the inner window.\n");
		if(nMain <= 0)
			die("nMain must be strictly positive.\n");
		if(nMain >= int(trialOrbitals.size()))
			die("nMain must be strictly less than the number of Wannier centers (%d).\n", int(trialOrbitals.size()));
	}
	if(innerWindow)
	{	if(!outerWindow) die("Inner window requires that an outer window be specified.\n");
		if(eInnerMin<eOuterMin || eInnerMax>eOuterMax)
			die("Inner window must lie entirely within the outer window.\n");
	}
	if(!outerWindow) //fixed bands
	{	int bStop = bStart+trialOrbitals.size();
		if(bStart<0 || bStop>e->eInfo.nBands)
			die("Index range [%d,%d) of participating bands incompatible with available bands [0,%d).\n", bStart, bStop, e->eInfo.nBands);
	}
	//Initialize minimizer:
	switch(localizationMeasure)
	{	case LM_FiniteDifference: wmin = std::make_shared<WannierMinimizerFD>(*e, *this); break;
		case LM_RealSpace: wmin = std::make_shared<WannierMinimizerRS>(*e, *this); break;
	}
	wmin->initTransformDependent();
	Citations::add("Maximally-localized Wannier functions",
		"N. Marzari and D. Vanderbilt, Phys. Rev. B 56, 12847 (1997)");
}

void Wannier::saveMLWF()
{	wmin->saveMLWF();
}

string Wannier::getFilename(FilenameType fnType, string varName, int* spin) const
{	string fname;
	switch(fnType)
	{	case FilenameInit: fname = initFilename; break;
		case FilenameDump: fname = dumpFilename; break;
	}
	string spinSuffix;
	if(spin && e->eInfo.spinType==SpinZ)
		spinSuffix = (*spin)==0 ? "Up" : "Dn";
	fname.replace(fname.find("$VAR"), 4, varName+spinSuffix);
	return fname;
}
