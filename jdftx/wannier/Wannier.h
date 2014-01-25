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

#ifndef JDFTX_ELECTRONIC_WANNIER_H
#define JDFTX_ELECTRONIC_WANNIER_H

#include <core/MinimizeParams.h>
#include <electronic/DOS.h>
#include <memory>

//! Compute Maximally-Localized Wannier Functions
class Wannier
{
public:
	Wannier();
	void setup(const Everything& everything);
	
	struct Center
	{	vector3<> r; //!< guess for center of localized wannier function
		double a; //!< exponential decay length of nodeless hydrogenic orbital of current l
		DOS::Weight::OrbitalDesc orbitalDesc; //!< orbital code
	};
	std::vector<Center> centers; //!< group of centers
	
	int bStart; //index of lowest band included in Wannier determination (used only when no energy windows)
	double eOuterMin, eOuterMax; //!< outer energy window (outside which bands do not contribute)
	double eInnerMin, eInnerMax; //!< inner energy window (within which all bands used)
	bool outerWindow, innerWindow; //!< denotes which windows are available
	
	bool saveWfns; //!< whether to write wavefunctions
	string initFilename, dumpFilename; //!< filename patterns for input and output
	
	void saveMLWF(); //!< Output the Maximally-Localized Wannier Functions from current wavefunctions
	
	//! Get filename for varName, based on initFilename if init=true and on dumpFilename otherwise.
	//! Optionally include Up/Dn suffix if spin is non-null and calculation is polarized
	string getFilename(bool init, string varName, int* spin=0) const; 
	
private:
	const Everything* e;
	MinimizeParams minParams;
	std::shared_ptr<class WannierMinimizer> wmin; //!< opaque struct to minimizer
	friend class WannierMinimizer;
	friend class CommandWannierMinimize;
};

#endif // JDFTX_ELECTRONIC_WANNIER_H
