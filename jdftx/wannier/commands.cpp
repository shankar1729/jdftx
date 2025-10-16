/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <commands/minimize.h>
#include <wannier/Wannier.h>
#include <core/LatticeUtils.h>

enum WannierMember
{	WM_addAtomicOrbitals,
	WM_pinAtomicOrbitals,
	WM_ignoreSemiCore,
	WM_localizationMeasure,
	WM_bStart,
	WM_outerWindow,
	WM_innerWindow,
	WM_projectionThresholds,
	WM_frozenCenters,
	WM_saveWfns,
	WM_saveWfnsRealSpace,
	WM_saveMomenta,
	WM_saveSpin,
	WM_saveR,
	WM_saveRP,
	WM_saveZ,
	WM_slabWeight,
	WM_loadRotations,
	WM_eigsOverride,
	WM_numericalOrbitals,
	WM_numericalOrbitalsOffset,
	WM_phononSup,
	WM_rSmooth,
	WM_spinMode,
	WM_polar,
	WM_delim
};

EnumStringMap<WannierMember> wannierMemberMap
(	WM_addAtomicOrbitals, "addAtomicOrbitals",
	WM_pinAtomicOrbitals, "pinAtomicOrbitals",
	WM_ignoreSemiCore, "ignoreSemiCore",
	WM_localizationMeasure, "localizationMeasure",
	WM_bStart, "bStart",
	WM_outerWindow, "outerWindow",
	WM_innerWindow, "innerWindow",
	WM_projectionThresholds, "projectionThresholds",
	WM_frozenCenters, "frozenCenters",
	WM_saveWfns, "saveWfns",
	WM_saveWfnsRealSpace, "saveWfnsRealSpace",
	WM_saveMomenta, "saveMomenta",
	WM_saveSpin, "saveSpin",
	WM_saveR, "saveR",
	WM_saveRP, "saveRP",
	WM_saveZ, "saveZ",
	WM_slabWeight, "slabWeight",
	WM_loadRotations, "loadRotations",
	WM_eigsOverride, "eigsOverride",
	WM_numericalOrbitals, "numericalOrbitals",
	WM_numericalOrbitalsOffset, "numericalOrbitalsOffset",
	WM_phononSup, "phononSupercell",
	WM_rSmooth, "rSmooth",
	WM_spinMode, "spinMode",
	WM_polar, "polar"
);

EnumStringMap<Wannier::LocalizationMeasure> localizationMeasureMap
(	Wannier::LM_FiniteDifference, "FiniteDifference",
	Wannier::LM_RealSpace, "RealSpace"
);

EnumStringMap<Wannier::SpinMode> spinModeMap
(	Wannier::SpinUp, "Up",
	Wannier::SpinDn, "Dn",
	Wannier::SpinAll, "All"
);


struct CommandWannier : public Command
{
	CommandWannier() : Command("wannier", "wannier")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"%Control calculation and output of maximally-localized Wannier functions \\cite MLWF.\n"
			"The possible <key>'s and their corresponding arguments are:\n"
			"\n+ addAtomicOrbitals yes|no\n\n"
			"   Whether to automatically add atomic orbitals for each atom.\n"
			"   Typically, this should be used with no wannier-center commands,\n"
			"   and is recommended when using projection based windows.\n"
			"   Default: no.\n"
			"\n+ pinAtomicOrbitals yes|no\n\n"
			"   Whether to pin the atomic orbitals during minimization,\n"
			"   as in wannier-center-pinned. Default: no.\n"
			"\n+ ignoreSemiCore yes|no\n\n"
			"   Whether to drop inner orbitals of each angular momentum when adding\n"
			"   atomic orbitals automatically (no effect if addAtomicOrbitals = no).\n"
			"   Default: yes.\n"
			"\n+ localizationMeasure FiniteDifference | RealSpace\n\n"
			"   Controls how the localization of the %Wannier functions is calculated.\n"
			"   The finite-difference reciprocal space measure of Marzari and Vanderbilt\n"
			"   (default) is inexpensive, but its error scales as Nkpoints^(-2./3).\n"
			"   The real space version is slower but its error scales as exp(-Nkpoints^(1/3)),\n"
			"   and is preferable for quantitative applications. Note that the real-space version\n"
			"   is not translationally invariant and wraps on a superlattice Wigner-Seitz cell\n"
			"   centered at the origin.\n"
			"\n+ bStart <band>\n\n"
			"   For fixed band calculations, 0-based index of lowest band used.\n"
			"   If atomic orbitals have been added with semi-core ignored,\n"
			"   this index excludes the ignored semi-core orbitals.\n"
			"   The number of bands equals the number of wannier-centers specified.\n"
			"   Default: 0. Use in insulator calculations to ignore semi-core orbitals.\n"
			"\n+ outerWindow <eMin> <eMax>\n\n"
			"   Energy window within which bands contribute to wannier functions \\cite MLWFmetal.\n"
			"   bStart is ignored if outerWindow is specified.\n"
			"\n+ innerWindow <eMin> <eMax>\n\n"
			"   Inner energy window within which bands are used exactly \\cite MLWFmetal.\n"
			"   Requires outerWindow, and innerWindow must be its subset.\n"
			"\n+ projectionThresholds <pOuter> <pInner>\n\n"
			"   Bands with higher projection than <pOuter> will be included in outer window.\n"
			"   Bands with higher projection than <pInner> will be included in inner window.\n"
			"   <pOuter> must be smaller than <pInner>. If energy windows are also specified,\n"
			"   the effective windows are unions of the energy and projection-based selections.\n"
			"\n+ frozenCenters <nFrozen> <filename>\n\n"
			"   Include frozen %Wannier centers imported as unitary rotations from <filename>\n"
			"   (.mlwfU output from another wannier run on the same jdftx state), and force\n"
			"   the current centers to be orthogonal to the imported ones. Note that all output\n"
			"   of this run will include the frozen as well as new centers for convenience.\n"
			"   \n"
			"   Requires innerWindow if outerWindow is present, and the frozen centers\n"
			"   must lie entirely within the inner window or fixed band range. (Effectively,\n" 
			"   the outer window used in generating the frozen centers must be a subset of\n"
			"   the present inner window / fixed band range.)\n"
			"\n+ saveWfns yes|no\n\n"
			"   Whether to write supercell wavefunctions in a spherical reciprocal-space basis.\n"
			"   Default: no.\n"
			"\n+ saveWfnsRealSpace yes|no\n\n"
			"   Whether to write supercell wavefunctions band-by-band in real space (can be enormous).\n"
			"   Default: no.\n"
			"\n+ saveMomenta yes|no\n\n"
			"   Whether to write momentum matrix elements in the same format as Hamiltonian.\n"
			"   The output is real and antisymmetric, i.e. it is [r,H] instead of p.\n"
			"   Default: no.\n"
			"\n+ saveSpin yes|no\n\n"
			"   Whether to write spin matrix elements (for non-collinear calculations only).\n"
			"   Default: no.\n"
			"\n+ saveR yes|no\n\n"
			"   Whether to write R matrix elements. Note that the position R here is only\n"
			"   the short-ranged part; longer-ranged contributions from the derivative\n"
			"   of the Wannier rotations need to be added after Wannier interpolation.\n"
			"   Default: no.\n"
			"\n+ saveRP yes|no\n\n"
			"   Whether to write R*P matrix elements for calculation of angular momentum and\n"
			"   electric quadrupole matrix elements in post-processing. Note that the position\n"
			"   R here is only the short-ranged part; longer-ranged contributions from\n"
			"   derivative of Wannier rotations need to be added after Wannier interpolation.\n"
			"   The output drops a factor of (-i) to make it real when possible, i.e.\n"
			"   it actually corresponds to R*[R,H] instead of R*P.\n"
			"   Default: no.\n"
			"\n+ saveZ <Vfilename> <Emag>\n\n"
			"   If specified, output matrix elements of z for perturbative electric field in post processing.\n"
			"   <Vfilename> is Vscloc output from a calculation with an applied electric field, and it\n"
			"   should be in the dump-name format contain $VAR (as in command fix-electron-potential).\n"
			"   <Emag> is the magnitude of the difference in electric field (in Eh/a0) between that calculation\n"
			"   and the present case. This combination is used to compute the matrix element accounting\n"
			"   for static local field effects due to the response of the system to the electric field.\n"
			"   Alternately, specify <Vfilename>=Ramp and <Emag>=0 to use a ramp potential (no local field).\n"
			"   Default: no.\n"
			"\n+ slabWeight <z0> <zH> <zSigma>\n\n"
			"   If specified, output the Wannier matrix elements of a slab weight function\n"
			"   centered at z0 (lattice coordinates) with half-width zH (lattice coordinates)\n"
			"   and Gaussian/erfc smoothness zSigma in bohrs. (Default: none)\n"
			"\n+ loadRotations yes|no\n\n"
			"   Whether to load rotations (.mlwU) from a previous %Wannier run.\n"
			"   Default: no.\n"
			"\n+ eigsOverride <filename>\n\n"
			"   Optionally read an alternate eigenvalues file to over-ride those from the total\n"
			"   energy calculation. Useful for generating Wannier Hamiltonians using eigenvalues\n"
			"   from another DFT or a many-body perturbation theory code.\n"
			"\n+ numericalOrbitals <filename>\n\n"
			"   Load numerical orbitals from <filename> with basis described in <filename>.header\n"
			"   that can then be used as trial orbitals. The reciprocal space wavefunction output\n"
			"   of wannier is precisely in this format, so that supercell %Wannier calculations can\n"
			"   be initialized from previously calculated bulk / unit cell %Wannier functions.\n"
			"\n+ numericalOrbitalsOffset <x0> <x1> <x2>\n\n"
			"   Origin of the numerical orbitals in lattice coordinates of the input.\n"
			"   The default [ .5 .5 .5 ] (cell center) is appropriate for output from wannier.\n"
			"\n+ phononSupercell <N0> <N1> <N2>\n\n"
			"   If specified, wannier will read in output (phononHsub and phononBasis) of the\n"
			"   phonon code with this supercell and output Wannierized electron-phonon matrix\n"
			"   elements (mlwfHePh). This file will contain matrices for pairs of cells in the\n"
			"   order specified in the mlwfCellMapSqPh output file, with an outer loop over\n"
			"   the nuclear displacement modes in phononBasis.\n"
			"\n+ rSmooth <rSmooth>\n\n"
			"   Width in bohrs of the supercell boundary region over which matrix elements are smoothed.\n"
			"   If phononSupercell is specified to process phonon quantities, the rSmooth specified here\n"
			"   must exactly match the value specified in the calculation in command phonon.\n"
			"\n+ spinMode" + spinModeMap.optionList() + "\n\n"
			"   If Up or Dn, only generate Wannier functions for that spin channel, allowing\n"
			"   different input files for each channel (independent centers, windows etc.).\n"
			"   Default: All, generate for all spin channels (only valid option for spintype != z-spin).\n"
			"\n+ polar" + boolMap.optionList() + "\n\n"
			"   Whether to include polar contribution subtractions in electron-phonon matrix elements.\n"
			"   Requires files 'totalE.Zeff' containing Born effective charges and 'totalE.epsInf'\n"
			"   containing optical dielectric tensor computed externally.\n"
			"   Default: false.";
		
		require("spintype");
		require("coulomb-interaction");
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		while(true)
		{	WannierMember key; pl.get(key, WM_delim, wannierMemberMap, "key");
			if(key==WM_delim) break;
			switch(key)
			{	case WM_addAtomicOrbitals:
					pl.get(wannier.addAtomicOrbitals, false,  boolMap, "addAtomicOrbitals", true);
					break;
				case WM_pinAtomicOrbitals:
					pl.get(wannier.pinAtomicOrbitals, false,  boolMap, "pinAtomicOrbitals", true);
					break;
				case WM_ignoreSemiCore:
					pl.get(wannier.ignoreSemiCore, true,  boolMap, "ignoreSemiCore", true);
					break;
				case WM_localizationMeasure:
					pl.get(wannier.localizationMeasure, Wannier::LM_FiniteDifference,  localizationMeasureMap, "localizationMeasure", true);
					break;
				case WM_bStart:
					pl.get(wannier.bStart, 0, "band", true);
					break;
				case WM_outerWindow:
					pl.get(wannier.eOuterMin, 0., "eMin", true);
					pl.get(wannier.eOuterMax, 0., "eMax", true);
					wannier.outerWindow = true;
					break;
				case WM_innerWindow:
					pl.get(wannier.eInnerMin, 0., "eMin", true);
					pl.get(wannier.eInnerMax, 0., "eMax", true);
					wannier.innerWindow = true;
					break;
				case WM_projectionThresholds:
					pl.get(wannier.projectionOuter, 0., "pOuter", true);
					pl.get(wannier.projectionInner, 1., "pInner", true);
					wannier.useProjectionThresholds = true;
					break;
				case WM_frozenCenters:
					pl.get(wannier.nFrozen, 0, "nFrozen", true);
					pl.get(wannier.frozenUfilename, string(), "filename", true);
					if(wannier.nFrozen <= 0) throw string("nFrozen must be > 0");
					break;
				case WM_saveWfns:
					pl.get(wannier.saveWfns, false, boolMap, "saveWfns", true);
					break;
				case WM_saveWfnsRealSpace:
					pl.get(wannier.saveWfnsRealSpace, false, boolMap, "saveWfnsRealSpace", true);
					break;
				case WM_saveMomenta:
					pl.get(wannier.saveMomenta, false, boolMap, "saveMomenta", true);
					break;
				case WM_saveSpin:
					pl.get(wannier.saveSpin, false, boolMap, "saveSpin", true);
					if(wannier.saveSpin and not e.eInfo.isNoncollinear())
						throw string("saveSpin requires noncollinear spin mode");
					break;
				case WM_saveR:
					pl.get(wannier.saveR, false, boolMap, "saveR", true);
					break;
				case WM_saveRP:
					pl.get(wannier.saveRP, false, boolMap, "saveRP", true);
					break;
				case WM_saveZ:
					pl.get(wannier.zVfilename, string(), "Vfilename", true);
					pl.get(wannier.zFieldMag, 0., "Emag", true);
					if(wannier.zVfilename == "Ramp")
					{	if(not e.coulombParams.isTruncated()[2])
							throw string("saveZ in ramp mode requires third lattice direction to be truncated");
					}
					else
					{	if(wannier.zVfilename.find("$VAR") == string::npos)
							throw string("<Vfilename> must contain $VAR");
						if(not wannier.zFieldMag)
							throw string("<Emag> must be non-zero when reading in another potential");
					}
					break;
				case WM_slabWeight:
					pl.get(wannier.z0, 0., "z0", true);
					pl.get(wannier.zH, 0., "zH", true);
					pl.get(wannier.zSigma, 0., "zSigma", true);
					break;
				case WM_loadRotations:
					pl.get(wannier.loadRotations, false, boolMap, "loadRotations", true);
					break;
				case WM_eigsOverride:
					pl.get(wannier.eigsFilename, string(), "filename", true);
					break;
				case WM_numericalOrbitals:
					pl.get(wannier.numericalOrbitalsFilename, string(), "filename", true);
					break;
				case WM_numericalOrbitalsOffset:
					pl.get(wannier.numericalOrbitalsOffset[0], 0., "x0", true);
					pl.get(wannier.numericalOrbitalsOffset[1], 0., "x1", true);
					pl.get(wannier.numericalOrbitalsOffset[2], 0., "x2", true);
					break;
				case WM_phononSup:
					pl.get(wannier.phononSup[0], 0, "N0", true);
					pl.get(wannier.phononSup[1], 0, "N1", true);
					pl.get(wannier.phononSup[2], 0, "N2", true);
					break;
				case WM_rSmooth:
					pl.get(wannier.rSmooth, 1., "rSmooth", true);
					if(wannier.rSmooth <= 0.) throw string("<rSmooth> must be positive");
					break;
				case WM_spinMode:
					pl.get(wannier.spinMode, Wannier::SpinAll,  spinModeMap, "spinMode", true);
					if(e.eInfo.spinType!=SpinZ && wannier.spinMode!=Wannier::SpinAll)
						throw string("<spinMode> must be All for spintype != z-spin");
					break;
				case WM_polar:
					pl.get(wannier.polar, false, boolMap, "polar", true);
					break;
				case WM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier& wannier = ((const WannierEverything&)e).wannier;
		logPrintf(" \\\n\taddAtomicOrbitals %s", boolMap.getString(wannier.addAtomicOrbitals));
		logPrintf(" \\\n\tpinAtomicOrbitals %s", boolMap.getString(wannier.pinAtomicOrbitals));
		logPrintf(" \\\n\tignoreSemiCore %s", boolMap.getString(wannier.ignoreSemiCore));
		logPrintf(" \\\n\tlocalizationMeasure %s", localizationMeasureMap.getString(wannier.localizationMeasure));
		logPrintf(" \\\n\tsaveWfns %s", boolMap.getString(wannier.saveWfns));
		logPrintf(" \\\n\tsaveWfnsRealSpace %s", boolMap.getString(wannier.saveWfnsRealSpace));
		logPrintf(" \\\n\tsaveMomenta %s", boolMap.getString(wannier.saveMomenta));
		logPrintf(" \\\n\tsaveSpin %s", boolMap.getString(wannier.saveSpin));
		logPrintf(" \\\n\tsaveR %s", boolMap.getString(wannier.saveR));
		logPrintf(" \\\n\tsaveRP %s", boolMap.getString(wannier.saveRP));
		if(wannier.zVfilename.length())
			logPrintf(" \\\n\tsaveZ %s %lg", wannier.zVfilename.c_str(), wannier.zFieldMag);
		if(wannier.zH)
			logPrintf(" \\\n\tslabWeight %lg %lg %lg", wannier.z0, wannier.zH, wannier.zSigma);
		logPrintf(" \\\n\tloadRotations %s", boolMap.getString(wannier.loadRotations));
		if(wannier.eigsFilename.length())
			logPrintf(" \\\n\teigsFilename %s", wannier.eigsFilename.c_str());
		
		if(wannier.outerWindow)
		{	logPrintf(" \\\n\touterWindow %lg %lg", wannier.eOuterMin, wannier.eOuterMax);
			if(wannier.innerWindow)
				logPrintf(" \\\n\tinnerWindow %lg %lg", wannier.eInnerMin, wannier.eInnerMax);
		}
		if(wannier.useProjectionThresholds)
			logPrintf(" \\\n\tprojectionThresholds %lg %lg", wannier.projectionOuter, wannier.projectionInner);
		if(not (wannier.outerWindow or wannier.useProjectionThresholds))
			logPrintf(" \\\n\tbStart %d", wannier.bStart);
		
		if(wannier.nFrozen)
			logPrintf(" \\\n\tfrozenCenters %d %s", wannier.nFrozen, wannier.frozenUfilename.c_str());
		if(wannier.numericalOrbitalsFilename.length())
		{	logPrintf(" \\\n\tnumericalOrbitals %s", wannier.numericalOrbitalsFilename.c_str());
			const vector3<>& offs = wannier.numericalOrbitalsOffset;
			logPrintf(" \\\n\tnumericalOrbitalsOffset %lg %lg %lg", offs[0], offs[1], offs[2]);
		}
		if(wannier.phononSup.length_squared())
			logPrintf(" \\\n\tphononSupercell %d %d %d", wannier.phononSup[0], wannier.phononSup[1], wannier.phononSup[2]);
		logPrintf(" \\\n\trSmooth %lg", wannier.rSmooth);
		logPrintf(" \\\n\tspinMode %s", spinModeMap.getString(wannier.spinMode));
		logPrintf(" \\\n\tpolar %s", boolMap.getString(wannier.polar));
	}
}
commandWannier;


struct CommandWannierCenter : public Command
{
	CommandWannierCenter(string suffix=string()) : Command("wannier-center" + suffix, "wannier")
	{
		format = "<aorb1> [<aorb2> ...]";
		comments =
			"Specify trial orbital for a wannier function as a linear combination of\n"
			"atomic orbitals <aorb*>. The syntax for each <aorb> is one of:\n"
			"\n"
			" +  <species> <atom> <orbDesc> [<coeff>=1.0]\n"
			" +  Gaussian <x0> <x1> <x2>  [<sigma>=1] [<orbDesc>=s] [<coeff>=1.0]\n"
			" +  Numerical <b>  [<x0> <x1> <x2>] [<coeff>=1.0]\n"
			"\n"
			"The first syntax selects an atomic orbital of the <atom>th atom\n"
			"(1-based index) of species named <species>.\n"
			"Orbital code <orbDesc> is as in command density-of-states,\n"
			"but only codes for individual orbitals eg. px, dxy can be used.\n"
			"(Codes for a set of orbitals eg. p, d are not allowed here.)\n"
			"\n"
			"The second syntax selects a Gaussian orbital of sigma = n*<sigma> bohrs,\n"
			"where n is the principal quantum number in <orbDesc> (default 1),\n"
			"and with angular quantum numbers as specified in <orbDesc> (default s).\n"
			"The orbital will be centered at <x0>,<x1>,<x2> in the coordinate system set by coords-type.\n"
			"\n"
			"The third syntax selects numerical orbital number <b> (0-based index)\n"
			"read from the file specified in command wannier, and optionally\n"
			"offsets it to the location specified by <x0>,<x1>,<x2>.\n"
			"\n"
			"When using multiple orbitals (of any type) in a linear combination,\n"
			"specify all default parameters explicitly before providing <coeff>.\n"
			"\n"
			"Specify this command once for each Wannier function.\n"
			"\n"
			"Examples:\n"
			"+ wannier-center Cu 1 dxy            #dxy orbital on first Cu atom.\n"
			"+ wannier-center Gaussian 0 0 0 1.5  #Gaussian orbital at origin with sigma=1.5.\n"
			"+ wannier-center Numerical 1         #Second numerical orbital without offset.\n"
			"+ wannier-center Cu 1 dxy 0.707  Numerical 1 0 0 0 0.707  #Linear combination of first and third example above.";
		
		allowMultiple = true;
		require("wannier-initial-state");
		require("wannier-dump-name");
		require("ion");
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		Wannier::TrialOrbital t;
		while(true)
		{	Wannier::AtomicOrbital ao;
			string key; pl.get(key, string(), "species");
			if(!key.size()) break;
			bool isGaussian = (key=="Gaussian");
			bool isNumerical = (key=="Numerical");
			if(isGaussian || isNumerical)
			{	if(isNumerical)
				{	pl.get(ao.numericalOrbIndex, -1, "b", true); //required orbital index
					if(ao.numericalOrbIndex<0) throw string("0-based numerical orbital index must be non-negative");
				}
				//Get position:
				pl.get(ao.r[0], 0., "x0", isGaussian); //position required for Gaussian, optional for numerical
				pl.get(ao.r[1], 0., "x1", isGaussian);
				pl.get(ao.r[2], 0., "x2", isGaussian);
				//Transform to lattice coordinates if necessary:
				if(e.iInfo.coordsType == CoordsCartesian)
					ao.r = inv(e.gInfo.R)*ao.r;
				if(isGaussian)
					pl.get(ao.sigma, 1., "sigma"); //optional sigma
			}
			else //must be atomic orbital
			{	for(int sp=0; sp<int(e.iInfo.species.size()); sp++)
					if(e.iInfo.species[sp]->name == key)
					{	ao.sp = sp;
						break;
					}
				if(ao.sp<0) //defaults to -1 in AtomicOrbital()
					throw string("<species> must match one of the atom-type names in the calculation, or 'Gaussian' or Numerical'");
				pl.get(ao.atom, 0, "atom", true);
				if(ao.atom <= 0) throw string("1-based atom index must be positive");
				ao.atom--; //store as zero-based index internally
				wannier.needAtomicOrbitals = true;
			}
			//Read orbital description if needed:
			if(!isNumerical) //i.e. Gaussian or atomic orbital
			{	string orbDesc;
				pl.get(orbDesc, string(), "orbDesc");
				if(orbDesc.length())
				{	ao.orbitalDesc.parse(orbDesc);
					if(ao.orbitalDesc.m > ao.orbitalDesc.l)
						throw string("Must specify a specific projection eg. px,py (not just p)");
				}
				else 
				{	if(isGaussian)
					{	//default s orbital
						ao.orbitalDesc.l = 0;
						ao.orbitalDesc.m = 0;
						ao.orbitalDesc.n = 0;
						ao.orbitalDesc.s = 0;
						ao.orbitalDesc.spinType = SpinNone;
					}
					else throw string("Must specify <orbDesc> explicitly for atomic orbitals.\n");
				}
				if((e.eInfo.spinType==SpinOrbit || e.eInfo.spinType==SpinVector) && ao.orbitalDesc.spinType==SpinNone)
					throw string("Must specify an explicit spin projection for noncollinear modes");
				if((e.eInfo.spinType==SpinNone || e.eInfo.spinType==SpinZ) && ao.orbitalDesc.spinType!=SpinNone)
					throw string("Orbital projections must not specify spin for collinear modes");
				if(ao.orbitalDesc.spinType==SpinOrbit && ao.atom<0)
					throw string("Relativistic (l,j,mj) orbital projections must be centered on an atom");
			}
			pl.get(ao.coeff, 1., "coeff");
			t.push_back(ao);
		}
		if(!t.size()) throw(string("Trial orbital for each center must contain at least one atomic orbital"));
		wannier.trialOrbitals.push_back(t);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier& wannier = ((const WannierEverything&)e).wannier;
		const Wannier::TrialOrbital& t = wannier.trialOrbitals[iRep];
		for(const Wannier::AtomicOrbital& ao: t)
		{	if(t.size()>1) logPrintf(" \\\n\t");
			bool isNumerical = (ao.numericalOrbIndex>=0);
			bool isGaussian = (ao.sigma>0.);
			if(isGaussian || isNumerical)
			{	if(isNumerical)
					logPrintf("Numerical %d", ao.numericalOrbIndex);
				else
					logPrintf("Gaussian");
				vector3<> r = ao.r;
				if(e.iInfo.coordsType == CoordsCartesian)
					r = e.gInfo.R * r; //report cartesian positions
				logPrintf(" %lg %lg %lg", r[0], r[1], r[2]);
				if(isGaussian)
					logPrintf(" %lg", ao.sigma);
			}
			else //atomic
			{	logPrintf("%s %d", e.iInfo.species[ao.sp]->name.c_str(), ao.atom+1);
			}
			if(!isNumerical) //Gaussian or atomic:
			{	logPrintf(" %s", string(ao.orbitalDesc).c_str());
			}
			logPrintf("  %lg", ao.coeff);
		}
	}
}
commandWannierCenter;


struct CommandWannierCenterPinned : public CommandWannierCenter
{
	CommandWannierCenterPinned(string suffix=string()) : CommandWannierCenter("-pinned")
	{	comments =
			"Same as command wannier-center, except that the minimizer tries to\n"
			"center this orbital at the specified guess positions (weighted mean\n"
			"of the orbital centers, if using a  linear combination). This can help\n"
			"improve the minimizer conditioning if the centers are known by symmetry.";
	}
	
	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		CommandWannierCenter::process(pl, e);
		Wannier::TrialOrbital& t = wannier.trialOrbitals.back();
		t.pinned = true;
	}
}
commandWannierCenterPinned;


struct CommandWannierMinimize : public CommandMinimize
{	CommandWannierMinimize() : CommandMinimize("wannier", "wannier") {}
    MinimizeParams& target(Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		return wannier.minParams;
	}
    void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		wannier.minParams.energyDiffThreshold = 1e-8; //override default value (0.) in MinimizeParams.h
		wannier.minParams.dirUpdateScheme = MinimizeParams::LBFGS;
		CommandMinimize::process(pl, e);
	}
}
commandWannierMinimize;


struct CommandWannierFilenames : public Command
{	virtual string& getTarget(Everything&)=0; //derived class determines where to save the file
	
	CommandWannierFilenames(string cmdSuffix, string comment) : Command("wannier-"+cmdSuffix, "wannier")
	{	format = "<format>";
		comments = 
			"Control the filename pattern for wannier " + comment + ", where <format> must contain\n"
			"a single instance of $VAR, which will be substituted by the name of each variable.";
	}

	void process(ParamList& pl, Everything& e)
	{	string& target = getTarget(e);
		pl.get(target, string("$STAMP.$VAR"), "format");
		if(target.find("$VAR")==string::npos)
			throw "<format> = " + target + " doesn't contain the pattern $VAR";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", getTarget(e).c_str());
	}
};

struct CommandWannierInitialState : public CommandWannierFilenames
{	CommandWannierInitialState() : CommandWannierFilenames("initial-state", "state initialization") {}
	string& getTarget(Everything& e) { return ((WannierEverything&)e).wannier.initFilename; }
}
commandWannierInitialState;

struct CommandWannierDumpName : public CommandWannierFilenames
{	CommandWannierDumpName() : CommandWannierFilenames("dump-name", "output") {}
	string& getTarget(Everything& e) { return ((WannierEverything&)e).wannier.dumpFilename; }
}
commandWannierDumpName;




struct CommandDefectSupercell : public Command
{
	CommandDefectSupercell() : Command("defect-supercell", "wannier")
	{
		format = "<name> <sup0> <sup1> <sup2> <inFile> <center0> <center1> <center2>  <q> [<alignWidth>=5. <alignSmooth>=1.]";
		comments =
			"Export wannierized matrix elements to mlwfHd_<name> using a <sup0> x <sup1> x <sup2>\n"
			"supercell (must evenly divide the corresponding k-point sampling in each direction)\n"
			"based on a defect calculation specified by <inFile>. The defect calculation need not\n"
			"be on the same supercell as the output, and this supercell is determined automatically.\n"
			"\n"
			"The defect is assumed to be centered on location <center0> <center1> <center2>\n"
			"in the coordinate system specified by command coords-type in <inFile>. (For coords-type\n"
			"lattice, these coordinates are fractional with respect to the defect calculation.)\n"
			"The defect charge <q> (specified as excess electrons, similar to command charged-defect)\n"
			"is used to correct for long-ranged potential behavior as appropriate.\n"
			"\n"
			"Note that the reference supercell calculation must export Vscloc.\n"
			"Only norm-conserving pseudopotentials are supported for defect matrix-element evaluation.\n"
			"\n"
			"Optional <alignWidth> and <alignSmooth> control the width and smoothness (in bohrs) of\n"
			"the region used for electrostatic potential alignment between the supercell and unit cell.\n"
			"The region is specified relative to the boundary of the supercell Wigner-Seitz cell.\n"
			"Use <alignWidth> <= 0. to skip the alignment potential altogether.\n"
			"\n"
			"This command may be specified multiple times to calculate elements for several defects.";
		
		require("coords-type");
		require("latt-scale");
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		DefectSupercell ds;
		pl.get(ds.name, string(), "name", true);
		for(int dir=0; dir<3; dir++) pl.get(ds.supOut[dir], 0, "sup"+string(1,"012"[dir]), true);
		pl.get(ds.inFile, string(), "inFile", true);
		for(int dir=0; dir<3; dir++) pl.get(ds.xCenter[dir], 0., "center"+string(1,"012"[dir]), true);
		pl.get(ds.q, 0., "q", true);
		pl.get(ds.alignWidth, 5., "alignWidth");
		pl.get(ds.alignSmooth, 1., "alignSmooth");
		wannier.defects.push_back(ds);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier& wannier = ((const WannierEverything&)e).wannier;
		const DefectSupercell& ds = wannier.defects[iRep];
		logPrintf("%s  %d %d %d  %s  %lg %lg %lg  %lg  %lg %lg", ds.name.c_str(),
			ds.supOut[0], ds.supOut[1], ds.supOut[2], ds.inFile.c_str(),
			ds.xCenter[0], ds.xCenter[1], ds.xCenter[2], ds.q,
			ds.alignWidth, ds.alignSmooth);
	}
}
commandDefectSupercell;
