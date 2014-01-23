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
#include <wannier/WannierMinimizer.h>

void Wannier::setup(const Everything& everything)
{	e = &everything;
	//Initialize minimization parameters:
	minParams.fpLog = globalLog;
	minParams.linePrefix = "WannierMinimize: ";
	minParams.energyLabel = "rVariance";
	//Initialize wminuator:
	wmin = std::make_shared<WannierMinimizer>(*e, *this);
	Citations::add("Maximally-localized Wannier functions",
		"N. Marzari and D. Vanderbilt, Phys. Rev. B 56, 12847 (1997)");
}

void Wannier::saveMLWF()
{	for(auto g: group)
		wmin->saveMLWF(g);
}

string Wannier::getFilename(bool init, string varName) const
{	string fname = init ? initFilename : dumpFilename;
	fname.replace(fname.find("$VAR"), 4, varName);
	return fname;
}
