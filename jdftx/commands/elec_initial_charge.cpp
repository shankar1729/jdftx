/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
#include <limits>

struct CommandElecInitialCharge : public Command
{
    CommandElecInitialCharge() : Command("elec-initial-charge")
	{
		format = "<QNet> | <QupNet> <QdnNet>";
		comments =
			"Initialize a system with <QNet> excess unpolarized electrons\n"
			"or with <QupNet> and <QdnNet> excess electrons in each spin channel.\n"
			"Both versions are valid irrespective of whether system is spin-polarized.";

		forbid("target-mu");
	}

	void process(ParamList& pl, Everything& e)
	{	double qTemp;
		e.eInfo.qNet.clear();
		//First q (required):
		pl.get(qTemp, 0., "QNet", true);
		e.eInfo.qNet.push_back(qTemp);
		//Second q (optional):
		pl.get(qTemp, std::numeric_limits<double>::quiet_NaN(), "QNet");
		if(!isnan(qTemp)) e.eInfo.qNet.push_back(qTemp);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%le", e.eInfo.qNet[0]);
		if(e.eInfo.qNet.size()==2)
			logPrintf(" %le", e.eInfo.qNet[1]);
	}
}
CommandElecInitialCharge;
