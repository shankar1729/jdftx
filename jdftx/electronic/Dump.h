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
#include <electronic/matrix.h>
#include <core/ScalarField.h>
#include <set>
#include <memory>

//! Dump frequency options:
enum DumpFrequency
{	DumpFreq_End, //!< At the end of calculation
	DumpFreq_Electronic, //!< Every (few) electronic step(s)
	DumpFreq_Fluid, //!< Every (few) fluid step(s)
	DumpFreq_Ionic, //!< Every (few) ionic step(s)
	DumpFreq_Lattice, //!< Every (few) lattice step(s)
	DumpFreq_Gummel, //!< Every (few) gummel step(s)
	DumpFreq_Delim  //special value used as a delimiter during command processing
};

//! Dump variable selection options:
enum DumpVariable { DumpAll, DumpNone, DumpState, //All, none or only those required to restart calculation
	DumpIonicPositions, DumpForces, DumpLattice, DumpIonicDensity, //Ionic positions, Forces, Lattice vectors, Nuclear charge density
	DumpElecDensity, DumpCoreDensity, DumpKEdensity, DumpFluidDensity, // electronic valence, core and KE densities, fluid densities
	DumpDvac, DumpDfluid, DumpDtot, //electrostatic potential of explicit system, fluid system, total
	DumpVcavity, DumpVfluidTot, //cavity potential of fluid, net electron potential due to fluid (electrostatic+cavity)
	DumpVlocps, DumpVscloc, DumpBandEigs, DumpEigStats, DumpFillings, DumpRhoAtom,
	DumpEcomponents, DumpExcCompare,
 	DumpBoundCharge, DumpSolvationRadii, DumpQMC, DumpOcean, DumpRealSpaceWfns, DumpFluidDebug, DumpSlabEpsilon, DumpChargedDefect,
	DumpOptVext, DumpDOS, DumpPolarizability, DumpElectronScattering, DumpSIC, DumpDipole, DumpStress, DumpExcitations,
	DumpMomenta, DumpSymmetries, DumpKpoints, DumpGvectors, DumpOrbitalDep, DumpXCanalysis, DumpEresolvedDensity, DumpFermiDensity,
	DumpDelim, //special value used as a delimiter during command processing
};


//! Stores the list of what to dump and when, and implements the functions to do so
class Dump : public std::set<std::pair<DumpFrequency,DumpVariable> >
{
public:
	void setup(const Everything&);
	
	//! Dump all variables that should be dumped at the given DumpFrequency type
	//!(End/Gummel/Ionic) at the iter'th iteration at that frequency
	void operator()(DumpFrequency freq, int iter);
	
	//! Get the dump filename corresponding to a particular variable name
	string getFilename(string varName) const;
	
	//! Check whether to dump at given frequency and iteration:
	bool checkInterval(DumpFrequency freq, int iter) const;
	
	std::shared_ptr<struct DOS> dos; //!< density-of-states calculator
	std::shared_ptr<struct Polarizability> polarizability; //!< electronic polarizability calculator
	std::shared_ptr<struct ElectronScattering> electronScattering; //!< electron-electron scattering calculator
	std::vector< std::pair<double,double> > densityErange; //!< energy ranges for energy-resolved density output
	std::vector<double> fermiDensityLevels; //!< energies at which to evaluate fermi-dirac derivative
	std::shared_ptr<struct SlabEpsilon> slabEpsilon; //!< slab dielectric function calculator
	std::shared_ptr<struct ChargedDefect> chargedDefect; //!< charged defect correction calculator
private:
	const Everything* e;
	string format; //!< Filename format containing $VAR, $STAMP, $FREQ etc.
	string stamp; //!< timestamp for current dump
	std::map<DumpFrequency,int> interval; //!< for each frequency, dump every interval times
	friend class CommandDump;
	friend class CommandDumpName;
	friend class CommandDumpInterval;
	void dumpQMC(); //!< QMC export implemented in DumpQMC.cpp
	void dumpOcean(); //!< BSE code export implemented in DumpOcean.cpp
	void dumpRsol(ScalarField nbound, string fname);
};

namespace Moments{void rn_pow_x(int i, vector3<> r, int dir, matrix3<> R, double moment, vector3<> r0, double* rx);}

#endif // JDFTX_ELECTRONIC_DUMP_H
