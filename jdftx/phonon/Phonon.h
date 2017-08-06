/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_PHONON_PHONON_H
#define JDFTX_PHONON_PHONON_H

#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/LatticeUtils.h>

//! @addtogroup Output
//! @{
//! @file Phonon.h Phonon class and related utilities

//! Block rotation matrix: used for symmetrization of electronic states
class BlockRotationMatrix
{	int nBlocks, blockSize;
	std::vector<int> colOffset; //!< column offset (in blocks) for the non-zero rotation block for each row block
	std::vector<matrix> rots; //!< rotation matrix for each row block
public:
	void init(int nBlocks, int blockSize); //!< initializes all blocks to zero
	void set(int rowBlock, int colBlock, const matrix& rot); //!< set one block
	void allReduce(); //!< collect over MPI
	
	matrix transform(const matrix& in) const; //!< return rot * in * dagger(rot)
};

//!Add reference to class Phonon to Everything (for use with the parser)
struct PhononEverything: public Everything
{	class Phonon& phonon;
	PhononEverything(class Phonon& phonon);
};

//! Calculate phonon dispersion, free energies and electron-phonon matrix elements
class Phonon
{
public:
	std::vector<std::pair<string,string> > input; //!< input file contents
	
	vector3<int> sup; //!< phonon supercell 
	double dr; //!< perturbation amplitude in Cartesian coordinates
	double T; //!< temperature for free energy estimation
	double Fcut; //!< fillings cutoff for optimizing number of bands
	double rSmooth; //!< supercell boundary width over which matrix elements are smoothed
	bool dryRun; //!< whether this is a dry run (test setup only; skip calculation)
	
	int iPerturbation; //!< if >=0, only run one supercell calculation
	bool collectPerturbations; //!< if true, collect results of previously computed perturbations (skips supercell SCF/Minimize)
	
	Phonon();
	void setup(bool printDefaults); //!< setup unit cell and basis modes for perturbations
	void dump(); //!< main calculations (sequence of supercell calculations) as well as output
	
	PhononEverything e; //!< data for original unit cell
private:
	PhononEverything eSupTemplate; //!< uninitialized version of eSup, with various flags later used to create eSup for each mode
	std::shared_ptr<PhononEverything> eSup; //!< supercell data for current perturbation

	int nSpins, nSpinor; //!< number of explicit spins and spinor length
	int nBandsOpt; //!< optimized number of bands, accounting for Fcut
	int prodSup; //!< number of unit cells in supercell
	std::vector<SpaceGroupOp> symSup; //!< Space group of unperturbed supercell (with translations restricted to unit cell)
	std::vector< matrix3<> > symSupCart; //!< Cartesian symmetry rotation matrices for unperturbed supercell
	std::vector< std::vector<BlockRotationMatrix> > stateRot; //!< Unitary rotation of states involved in gamma-point Hsub for each supercell symmetry operation
	
	//!Basis for phonon modes (not reduced by symmetries):
	struct Mode
	{	int sp, at; //!< species and atom number (within first unit cell)
		vector3<> dir; //!< Cartesian unit vector for displacement
	};
	std::vector<Mode> modes;
	double E0; //!< energy of unperturbed supercell
	IonicGradient grad0; //!< forces of unperturbed supercell
	std::vector<IonicGradient> dgrad; //!< change in forces per unit displacement in each mode (force matrix)
	std::vector< std::vector<matrix> > dHsub; //!< change in electronic subspace Hamiltonian per unit displacement in each mode (precursor to electron-phonon matrix elements)
	
	//!Minimal basis of perturbations for supercell calculations:
	struct Perturbation : public Mode
	{	double weight; //!< weight of perturbation (adds up to nModes / nSymmetries)
	};
	std::vector<Perturbation> perturbations;
	
	//!Run supercell calculation for specified perturbation (using fnamePattern to load/restore required properties)
	void processPerturbation(const Perturbation& pert, string fnamePattern);
	
	//!Set unperturbed state of supercell from unit cell and retrieve unperturbed subspace Hamiltonian at supercell Gamma point (for all bands)
	std::vector<diagMatrix> setSupState();
	
	//!Calculate subspace Hamiltonian of perturbed supercell:
	std::vector<matrix> getPerturbedHsub();
	
	//!Mapping between unit cell and current supercell k-points (generated by processPerturbation and used by setSupState)
	struct StateMapEntry : public Supercell::KmeshTransform //contain source k-point rotation here
	{	int qSup; //!< state index for supercell
		int nqPrev; //!< number of previous unit cell k-points that point to this supercell
		vector3<> k; //!< unit cell k-vector
		std::shared_ptr<class ColumnBundleTransform> transform; //!< wavefunction map
	};
	std::vector<StateMapEntry> stateMap; //!< map from unit cell k-points to supercell k-points
	
	//Utilities for mapping 3-dimensional and flattened cell indices:
	vector3<int> getCell(int unit) const;
	int getUnit(const vector3<int>& cell) const;
	
	//! Enforce hermiticity of force matrix
	void forceMatrixDaggerSymmetrize(std::vector<matrix>& omegaSq, const std::map<vector3<int>,matrix>& cellMap) const;

	//! Enforce translational invariance sum rule on force matrix
	void forceMatrixEnforceSumRule(std::vector<matrix>& omegaSq, const std::map<vector3<int>,matrix>& cellMap, const std::vector<double>& invsqrtM) const;
};

//! @}
#endif //JDFTX_PHONON_PHONON_H
