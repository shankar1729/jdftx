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
#include <wannier/DefectSupercell.h>
#include <memory>

//! @addtogroup Output
//! @{
//! @file Wannier.h

//! Compute Maximally-Localized Wannier Functions
class Wannier
{
public:
	Wannier();
	void setup(const Everything& everything);
	
	struct AtomicOrbital
	{	vector3<> r; //!< guess for center of localized wannier function
		double sigma; //!< width parameter of Gaussian orbital of current l
		int sp, atom; //!< species code (<0 if not using a pseudopotential atomic orbital) and atom number
		int numericalOrbIndex; //!< index to a numerical orbital (<0 if not using a numerical orbital)
		DOS::Weight::OrbitalDesc orbitalDesc; //!< orbital code
		double coeff; //!< coefficient (prefactor) in contribution to trial orbital (1 if only using a single orbital)	
		AtomicOrbital() : sigma(0.), sp(-1), atom(-1), numericalOrbIndex(-1), coeff(1.) {}
	};
	struct TrialOrbital : public std::vector<AtomicOrbital>
	{	bool pinned;
		vector3<> xCenter; //center of trial orbital in lattice coordinates (used for pinning and improved translation invariance correction)
		TrialOrbital() : pinned(false) {}
	};
	std::vector<TrialOrbital> trialOrbitals; //!< group of centers
	bool needAtomicOrbitals;
	
	bool addAtomicOrbitals; //!< whether to automatically add atomic orbitals for all atoms
	bool pinAtomicOrbitals; //!< whether to pin centers of automatically-added atomic orbitals
	bool ignoreSemiCore; //!< whether to ignore semi-core orbitals when adding atomic orbitals
	
	enum LocalizationMeasure
	{	LM_FiniteDifference,
		LM_RealSpace
	}
	localizationMeasure;
	
	int bStart; //index of lowest band included in Wannier determination (used only when no energy windows)
	double eOuterMin, eOuterMax; //!< outer energy window (outside which bands do not contribute)
	double eInnerMin, eInnerMax; //!< inner energy window (within which all bands used)
	double projectionOuter; //!< projection threshold for inclusion in outer window
	double projectionInner; //!< projection threshold for inclusion in inner window
	bool outerWindow, innerWindow, useProjectionThresholds; //!< denotes which windows are available
	
	int nFrozen; string frozenUfilename; //!< number of frozen centers, and the filename to read their rotations from
	int nCenters; //!< total number of centers, those being optimized and frozen
	
	bool saveWfns; //!< whether to write wavefunctions
	bool saveWfnsRealSpace; //!< whether to output Wannier functions band-by-band in real-space
	bool saveMomenta; //!< whether to output momentum matrix elements
	bool saveSpin; //!< whether to output spin matrix elements (non-collinear only)
	bool saveR; //!< whether to output position matrix elements R (short-range Wannier-rotation-matched part)
	bool saveRP; //!< whether to output R*P matrix elements for Wannierized angular momentum / electric quadrupole evaluation
	
	string zVfilename; //!< filename for reading Vscloc with an applied electric field for z matrix element output
	double zFieldMag; //!< magnitude of electric field difference (Eh/a0) between current calculation and the specified Vscloc
	
	double z0, zH, zSigma; //!< center (lattice coords), half-width (lattice coords) and smoothness (bohrs) for slab-weight function

	bool loadRotations; //!< whether to load initial rotations from previous dump
	string initFilename, dumpFilename; //!< filename patterns for input and output
	string eigsFilename; //!< optional override for eigenvals file
	
	string numericalOrbitalsFilename; //!< filename for reading numerical orbitals
	vector3<> numericalOrbitalsOffset; //!< lattice coordinates of the origin in the input
	
	vector3<int> phononSup; //!< phonon supercell (process e-ph matrix elements on this supercell if non-zero)
	double rSmooth; //!< supercell boundary width over which matrix elements are smoothed
	
	enum SpinMode
	{	SpinUp,
		SpinDn,
		SpinAll
	}
	spinMode; //!< which spin(s) to generate Wannier functions for
	std::vector<int> iSpinArr; //!< set of spin indices corresponding to spinMode
	bool polar; //!< whether to subtract long-range contributions in e-ph matrix elements
	
	std::vector<DefectSupercell> defects; //!< List of defect supercells to compute Wannierized matrix elements for
	
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
	int nBandsSemiCore;
	friend class WannierMinimizer;
	friend class DefectSupercell;
	friend struct CommandWannierMinimize;
};

//!Version of Everything with Wannier added
struct WannierEverything : public Everything
{	Wannier wannier;
	void setup();
};

//! @}
#endif // JDFTX_WANNIER_WANNIER_H
