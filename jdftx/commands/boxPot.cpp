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

struct CommandBoxPot : public Command
{
	CommandBoxPot() : Command("box-potential")
	{
		format = "xmin xmax ymin ymax zmin zmax Vin Vout [<convolve_radius>=0.1]";
		comments =
			"Include an step-function shaped external potential (in hartrees) for the electrons";
	
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	
		double xmin, xmax, ymin, ymax, zmin, zmax;
		double Vin, Vout;
		double convolve_radius;
		
		pl.get(xmin, 0., "xmin", true);
		pl.get(xmax, 0., "xmax", true);
		pl.get(ymin, 0., "ymin", true);
		pl.get(ymax, 0., "ymax", true);
		pl.get(zmin, 0., "zmin", true);
		pl.get(zmax, 0., "zmax", true);
		
		if(xmin< 0 or xmax < xmin or ymin< 0 or ymax < ymin or xmin< 0 or zmax < zmin)
			die("max coordinates must be smaller than min coordinates with both being greater than 0");
		
		pl.get(Vin, 0., "zmax", true);
		pl.get(Vout, 0., "zmax", true);
		
		pl.get(convolve_radius, 0.2, "zmax", false);
		
		e.eVars.boxPot.push_back(BoxPotential(xmin, xmax, ymin, ymax, zmin, zmax, Vin, Vout, convolve_radius));

	}

	void printStatus(Everything& e, int iRep)
	{	auto bP = e.eVars.boxPot[iRep];
		logPrintf("%.5e %.5e %.5e %.5e %.5e %.5e    %.5e %.5e", bP.xmin, bP.xmax, bP.ymin, bP.ymax, bP.zmin, bP.zmax, bP.Vin, bP.Vout);
	}
}
commandBoxPot;
