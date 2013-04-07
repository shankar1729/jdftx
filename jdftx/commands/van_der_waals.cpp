/*-------------------------------------------------------------------
Copyright 2011 Deniz Gunceler

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
#include <electronic/VanDerWaals.h>

struct CommandVanDerWaals : public Command
{
	CommandVanDerWaals() : Command("van-der-waals")
	{
		format = "[<scaleOverride>=0]";
		comments = "Pair-potential corrections for the long range Van der Waals interaction.\n"
			"Implementation follows \"S. Grimme, J. Comput. Chem. 27: 1787–1799 (2006)\"\n"
			"Exchange-Correlation functionals supported with van-der-waals are gga-PBE, hyb-gga-xc-b3lyp and mgga-TPSS.\n"
			"Manually specify <scaleOverride> to use with other functionals";
	}

	void process(ParamList& pl, Everything& e)
	{	e.iInfo.vdWenable = true;
		pl.get(e.iInfo.vdWscale, 0., "scaleOverride");
	}

	void printStatus(Everything& e, int iRep)
	{	if(e.iInfo.vdWscale) logPrintf("%lg", e.iInfo.vdWscale);
	}
}
commandVanDerWaals;
