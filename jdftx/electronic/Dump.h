/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#ifndef JDFTX_ELECTRONIC_DUMP_H
#define JDFTX_ELECTRONIC_DUMP_H

#include <electronic/common.h>
#include <electronic/Wannier.h>
#include <set>
#include <memory>

//! Dump frequency options:
enum DumpFrequency
{	DumpFreq_End, //!< At the end of calculation
	DumpFreq_Electronic, //!< Every (few) electronic step(s)
	DumpFreq_Fluid, //!< Every (few) fluid step(s)
	DumpFreq_Ionic, //!< Every (few) ionic step(s)
	DumpFreq_Gummel, //!< Every (few) gummel step(s)
	DumpFreq_Delim  //special value used as a delimiter during command processing
};

//! Dump variable selection options:
enum DumpVariable { DumpAll, DumpNone, DumpState, //All, none or only those required to restart calculation
	DumpIonicPositions, DumpForces, DumpIonicDensity, //Ionic positions, Forces, Nuclear charge density
	DumpElecDensity, DumpCoreDensity, DumpFluidDensity, // electronic valence and core densities, fluid densities
	DumpDvac, DumpDfluid, DumpDtot, //electrostatic potential of explicit system, fluid system, total
	DumpVcavity, DumpVfluidTot, //cavity potential of fluid, net electron potential due to fluid (electrostatic+cavity)
	DumpVlocps, DumpVscloc, DumpHsubEvecs, DumpBandEigs,
	DumpEcomponents, DumpExcCompare,
	DumpBoundCharge, DumpQMC, DumpRealSpaceWfns, DumpFluidDebug,
	DumpSpinOrbit, DumpProjectors, DumpWannier, DumpOptVext, DumpDOS,
	DumpDelim //special value used as a delimiter during command processing
};


//! Stores the list of what to dump and when, and implements the functions to do so
class Dump : public std::set<std::pair<DumpFrequency,DumpVariable> >
{
public:
	void setup(const Everything&);
	
	//! Dump all variables that should be dumped at the given DumpFrequency type (End/Gummel/Ionic)
	void operator()(DumpFrequency freq);
	
	//! Get the dump filename corresponding to a particular variable name
	string getFilename(string varName) const;
	
	bool shouldDump(DumpFrequency freq, int iter) const; //!< whether dump of a particular frequency should happen at a given iteration

	Wannier wannier; //!< wannier function calculator
	std::shared_ptr<struct DOS> dos; //!< density-of-states calculator
private:
	const Everything* e;
	string format; //!< Filename format containing $VAR, $STAMP, $FREQ etc.
	string stamp; //!< timestamp for current dump
	std::map<DumpFrequency,int> interval; //!< for each frequency, dump every interval times
	friend class CommandDump;
	friend class CommandDumpName;
	friend class CommandDumpInterval;
	void dumpQMC(); //!< QMC export implemented in DumpQMC.cpp
};

#endif // JDFTX_ELECTRONIC_DUMP_H
