/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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

struct CommandCustomFilling : public Command
{
	int qnum, band; double filling;
	
	CommandCustomFilling() : Command("custom-filling")
	{
		format = "<qnum> <band> <filling>";
		comments = "Specify a custom filling for the input (quantum-number, band)\n"
					"Bands are indexed from HOMO, i.e. band=0 is HOMO, band=1 is the LUMO\n"
					"Fillings are in normalized units, as output by the dumpName.fillings file";
		allowMultiple = true;
		forbid("elec-initial-charge");
		forbid("initial-magnetic-moments");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(qnum, 0, "qnum", true);
		pl.get(band, 0, "band", true);
		pl.get(filling, 0., "filling", true);
		
		if((filling > (e.eInfo.spinType == SpinZ ? 1. : 2.)) or (filling < 0))
		{	ostringstream oss; oss << "Fillings must be between 0 and " << (e.eInfo.spinType == SpinZ ? 1. : 2.);
			throw oss.str();
		}
		
		e.eInfo.customFillings.push_back(std::make_tuple(qnum, band, filling));
		
	}

	void printStatus(Everything& e, int iRep)
	{	std::vector<std::tuple<int,int,double>>& customFillings = e.eInfo.customFillings;
		logPrintf("%i %i %.2e", std::get<0>(customFillings[iRep]), std::get<1>(customFillings[iRep]), std::get<2>(customFillings[iRep]));
	}
}
commandCustomFilling;