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

#ifndef JDFTX_FLUID_MOLECULE_H
#define JDFTX_FLUID_MOLECULE_H

#include <electronic/RadialFunction.h>
#include <core/vector3.h>
#include <core/string.h>
#include <core/ScalarField.h>
#include <vector>
#include <map>

//! Multi-site molecule model
struct Molecule
{	
	string name; //!< Molecule name
	
	struct Site
	{	string name; //!< site name
		double Rhs; //!< hard sphere radius
		int atomicNumber; //!< necessary for vdW parameters
		double Znuc, sigmaNuc; //!< magnitude of the nuclear charge (positive) and corresponding gaussian width
		double Zelec, aElec; //!< magnitude of electron charge (positive) and corresponding cuspless-exponential width
	        double Zsite; //!<site charge in electrons
		double deltaS; //!<G=0 correction potential due to electron-fluid and fluid-fluid charge kernel mismatch
		double sigmaElec, rcElec; //!< width and location of peak in electron charge distribution
		string elecFilename, elecFilenameG; //!< include electron charge density from real- or G- space radial ASCII file
		double alpha, aPol; //!< isotropic polarizability and corresponding cuspless-exponential width

		std::vector< vector3<> > positions; //!< Positions w.r.t molecular origin in the reference orientation

		Site(string name, int atomicNumber=0);
		~Site();
		void setup(const GridInfo& gInfo); //!< initialize the radial functions from the properties specified above
		explicit operator bool() const { return initialized; } //!< return whether site has been setup
		
		RadialFunctionG w0, w1, w2, w3, w1v, w2m; //!< Hard sphere weight functions
		RadialFunctionG elecKernel, chargeKernel, polKernel; //!< Electron density, net charge density and polarizability kernels for the sites
	private:
		bool initialized;
		void free();
	};
	std::vector< std::shared_ptr<Site> > sites;
	RadialFunctionG mfKernel; //!< Mean field interaction kernel (with minimum Coulomb self energy while preserving intermolecular interactions)

	Molecule(string name=string());
	~Molecule();
	void setup(const GridInfo& gInfo, double Rmf);
	explicit operator bool() const { return initialized; } //!< return whether site has been setup
	
	bool isMonoatomic() const; //!< whether it is a monoatomic molecule
	double getCharge() const; //!< total charge on molecule
	double checkCharge(); //!< check total charge on molecule and modify kernels if small non-neutral charge
	vector3<> getDipole() const; //!< total dipole moment on molecule
	double getVhs() const; //!< total exclusion volume
	double getAlphaTot() const; //!< total polarizability
	
	std::map<double,int> getBonds() const; //!< get the harmonic sum of radii for spheres in contact, with the multiplicities for each such pair
	
	void setModelMonoatomic(string name, double Q, double Rhs); //!< set to a simple monoatomic model (useful for debugging, not for actual solvation)
private:
	bool initialized;
};

#endif // JDFTX_FLUID_MOLECULE_H
