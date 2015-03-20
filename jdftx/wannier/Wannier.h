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

#ifndef JDFTX_WANNIER_WANNIER_H
#define JDFTX_WANNIER_WANNIER_H

#include <core/MinimizeParams.h>
#include <electronic/Everything.h>
#include <electronic/DOS.h>
#include <memory>

//! Compute Maximally-Localized Wannier Functions
class Wannier
{
public:
	Wannier();
	void setup(const Everything& everything);
	
	struct AtomicOrbital
	{	vector3<> r; //!< guess for center of localized wannier function
		double a; //!< exponential decay length of nodeless hydrogenic orbital of current l
		int sp, atom; //!< species code (<0 if not using a pseudopotential atomic orbital) and atom number (<0 if not using an orbital located on an atom)
		int numericalOrbIndex; //!< index to a numerical orbital (<0 if not using a numerical orbital)
		DOS::Weight::OrbitalDesc orbitalDesc; //!< orbital code
		double coeff; //!< coefficient (prefactor) in contribution to trial orbital (1 if only using a single orbital)
	};
	struct TrialOrbital : public std::vector<AtomicOrbital>
	{	bool pinned;
		vector3<> rCenter;
		TrialOrbital() : pinned(false) {}
	};
	std::vector<TrialOrbital> trialOrbitals; //!< group of centers
	bool needAtomicOrbitals;
	
	enum LocalizationMeasure
	{	LM_FiniteDifference,
		LM_RealSpace
	}
	localizationMeasure;
	bool precond; //whether to precondition minimize
	
	int bStart; //index of lowest band included in Wannier determination (used only when no energy windows)
	double eOuterMin, eOuterMax; //!< outer energy window (outside which bands do not contribute)
	double eInnerMin, eInnerMax; //!< inner energy window (within which all bands used)
	bool outerWindow, innerWindow; //!< denotes which windows are available
	
	int nFrozen; string frozenUfilename; //!< number of frozen centers, and the filename to read their rotations from
	int nCenters; //!< total number of centers, those being optimized and frozen
	
	bool saveWfns; //!< whether to write wavefunctions
	bool saveWfnsRealSpace; //!< whether to output Wannier functions band-by-band in real-space
	bool saveMomenta; //!< whether to output momentum matrix elements
	bool loadRotations; //!< whether to load initial rotations from previous dump
	string initFilename, dumpFilename; //!< filename patterns for input and output
	
	string numericalOrbitalsFilename; //!< filename for reading numerical orbitals
	vector3<> numericalOrbitalsOffset; //!< lattice coordinates of the origin in the input
	
	vector3<int> phononSup; //!< phonon supercell (process e-ph matrix elements on this supercell if non-zero)
	
	void saveMLWF(); //!< Output the Maximally-Localized Wannier Functions from current wavefunctions
	
	enum FilenameType
	{	FilenameInit,
		FilenameDump
	};
	//! Get filename for varName, based on initFilename, dumpFilename or numericalOrbitalsFilename depending on fnType
	//! Optionally include Up/Dn suffix if spin is non-null and calculation is polarized
	string getFilename(FilenameType fnType, string varName, int* spin=0) const; 
	
private:
	const Everything* e;
	MinimizeParams minParams;
	std::shared_ptr<class WannierMinimizer> wmin; //!< opaque struct to minimizer
	friend class WannierMinimizer;
	friend class CommandWannierMinimize;
};

//!Version of Everything with Wannier added
struct WannierEverything : public Everything
{	Wannier wannier;
	void setup();
};

#endif // JDFTX_WANNIER_WANNIER_H
