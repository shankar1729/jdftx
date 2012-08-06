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

#include <commands/command.h>
#include <electronic/Everything.h>
#include <electronic/InverseKohnSham.h>

struct CommandInvertKohnSham : public Command
{
	CommandInvertKohnSham() : Command("invertKohnSham")
	{
		format = "[<nonlocal>=yes] [<sigma>=0] [<chiGuessFilename>]";
		comments =
			"Solve inverse Kohn-Sham problem: for a given electron density\n"
			"specified using fix-electron-density, find the corresponding\n"
			"external potential (in addition to the pseudopotentials and\n"
			"Hartree+XC potential evaluated at the given electron density).\n"
			"Vexternal may be used to specify an initial guess.\n"
			"Control outer optimization using inverseKohnSham-minimize, and\n"
			"inner minimization using electronic-minimize (as usual). The\n"
			"result is dumped with variable name \"optVext\" or \"optVextUp\"\n"
			"and \"optVextDn\" depending on specified spin type.\n"
			"Option <nonlocal>="+boolMap.optionList()+" controls whether to\n"
			"include non-local parts of the pseudopotential (default yes).\n"
			"Option <sigma> specifies a bandwidth cutoff in the external\n"
			"potential of the form exp(-(1/2) sigma^2 G^2) (default: 0).\n"
			"Option <chiGuessFilename> specifies a preconditioner based on\n"
			"the response function of a similar electronic system. The pattern\n"
			"should contain $VAR which will be used to read wavefunctions,\n"
			"eigenvalues and fillings (these should include empty states).";
		
		require("fix-electron-density");
	}

	void process(ParamList& pl, Everything& e)
	{	e.cntrl.invertKS = true;
		e.dump.insert(std::make_pair(DumpFreq_End, DumpOptVext));
		//Optional parameters:
		pl.get(e.cntrl.invertKS_nonlocal, true, boolMap, "nonlocal");
		pl.get(e.cntrl.invertKS_sigma, 0., "sigma");
		string& fname = e.cntrl.invertKS_chiGuessFilename;
		pl.get(fname, string(), "chiGuessFilename");
		if(fname.length() && fname.find("$VAR")==string::npos)
			throw "<chiGuessFilename> = " + fname + " doesn't contain '$VAR'";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s %lg %s", boolMap.getString(e.cntrl.invertKS_nonlocal),
			e.cntrl.invertKS_sigma, e.cntrl.invertKS_chiGuessFilename.c_str());
	}
} commandInvertKohnSham;
