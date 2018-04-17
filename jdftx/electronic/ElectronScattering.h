/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_ELECTRONSCATTERING_H
#define JDFTX_ELECTRONIC_ELECTRONSCATTERING_H

#include <electronic/Basis.h>
#include <core/LatticeUtils.h>
#include <memory>

class ColumnBundle;
class diagMatrix;
class QuantumNumber;

//! @addtogroup Output
//! @{

//! Electron-electron scattering (ImSigma_ee) calculator
struct ElectronScattering
{
	double eta; //!< frequency resolution and half the imaginary part ascribed to probe frequency (set to eInfo.kT if unspecified)
	double Ecut; //!< energy cut-off for dielectric matrices (set to cntrl.Ecut if unspecified)
	double fCut; //!< threshold for considering states fully occupied / unoccupied (default: 1e-6)
	double omegaMax; //!< maximum energy transfer to account for and hence maximum frequency in dielectric grid (if zero, autodetermine from available eigenvalues)
	
	bool slabResponse; //!< whether to work in slab response output mode
	double EcutTransverse; //!< energy cutoff in directions transverse to slab normal (same as Ecut above if unspecified)
	
	ElectronScattering();
	void dump(const Everything& e); //!< compute and dump Im(Sigma_ee) for each eigenstate

private:
	const Everything* e;
	int nBands, nSpinor, nSpins, qCount;
	double Emin, Emax; //!< energy range that contributes to transitions less than omegaMax
	std::vector<ColumnBundle> C; //wavefunctions, made available on all processes
	std::vector<diagMatrix> E, F; //energies and fillings, available on all processes
	std::vector<std::vector<matrix>> VdagC; //pseudopotential projections
	std::shared_ptr<const Supercell> supercell; //contains transformations between full and reduced k-mesh
	std::shared_ptr<const PeriodicLookup< vector3<> > > plook; //O(1) lookup for finding k-points in mesh
	std::vector<QuantumNumber> qmesh; //reduced momentum-transfer mesh
	std::vector<Basis> basisChi; //polarizability bases for qmesh
	Basis basis; //common wavefunction  basis
	std::map< vector3<int>, std::shared_ptr<class ColumnBundleTransform> > transform; //k-mesh transformations
	std::map< vector3<int>, QuantumNumber > qnumMesh; //equivalent of eInfo.qnums for entire k-mesh
	std::vector<matrix> nAugRhoAtom; //augmentation pair densities per atomic density-matrix element for each species
	
	struct Event
	{	int i, j; //band indices
		double fWeight; //relevant fillings combination ( (fi - fj)/2 or (1-fi-fj) for chiMode = true or false )
		double Eji; //Ej - Ei
	};
	
	std::vector<Event> getEvents( //!< return list of events
		bool chiMode, //!< whether to calculate fWeight for chi (true) or ImSigma (false), and hence which events to select (occ-unocc or occ-occ + unocc-unocc)
		int iSpin, //!< spin index (0 or 1 if z-spin, 0 otherwise)
		size_t ik, //!< k-mesh index for first state
		size_t iq, //!< q-mesh index for momentum transfer
		size_t& jk, //!< set k-mesh index for second state
		matrix& nij //!< set pair densities for each event, one per column
	) const;
	
	ColumnBundle getWfns(size_t ik, int iSpin, const vector3<>& k, std::vector<matrix>* VdagCi=0) const; //get wavefunctions at an arbitrary point in k-mesh
	matrix coulombMatrix(size_t iq) const; //retrieve the Coulomb operator for a specific momentum transfer
	void nAugRhoAtomInit(size_t iq); //Initialize nAugRhoAtom for a specific momentum transfer
	void dumpSlabResponse(Everything& e, const diagMatrix& omegaGrid);
};

//! @}
#endif //JDFTX_ELECTRONIC_ELECTRONSCATTERING_H
