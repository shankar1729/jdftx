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

#ifndef JDFTX_ELECTRONIC_SYMMETRIES_H
#define JDFTX_ELECTRONIC_SYMMETRIES_H

#include <core/ScalarFieldArray.h>
#include <list>

class Everything;
class QuantumNumber;

//! @addtogroup ElectronicDFT
//! @{

enum SymmetryMode {SymmetriesNone, SymmetriesAutomatic, SymmetriesManual}; //!< symmetry modes

//! Symmetry detection and symmetrization of various quantities
class Symmetries
{
public:
	SymmetryMode mode; //!< Symmetry mode (none, automatic or manual)

	Symmetries();
	void setup(const Everything& everything); //!< Phase 1 of setup which computes/checks lattice+basis symmetries
	void setupMesh(); //!< Phase 2 of setup which computes / checks FFTbox and k-mesh dependent symmetries
	
	//! Reduce a k-point mesh (and remember its inversion symmetry property in kpointInvertList)
	std::vector<QuantumNumber> reduceKmesh(const std::vector<QuantumNumber>& qnums) const;
	
	void symmetrize(ScalarField&) const; //!< symmetrize a scalar field in real space
	void symmetrize(ScalarFieldTilde&) const; //!< symmetrize a scalar field in reciprocal space
	void symmetrize(complexScalarFieldTilde&) const; //!< symmetrize a complex scalar field in reciprocal space
	void symmetrize(ScalarFieldArray&) const; //!< symmetrize an array of scalar fields in real space representing spin density / potentials
	void symmetrize(ScalarFieldTildeArray&) const; //!< symmetrize an array of scalar fields in reciprocal space representing spin density / potentials
	void symmetrize(std::vector<complexScalarFieldTilde>&) const; //!< symmetrize an array of complex scalar fields in reciprocal space representing spin density / potentials
	void symmetrize(struct IonicGradient&) const; //!< symmetrize forces
	void symmetrizeSpherical(matrix&, const class SpeciesInfo* specie) const; //!< symmetrize matrices in Ylm basis per atom of species sp (accounting for atom maps)
	const std::vector<SpaceGroupOp>& getMatrices() const; //!< directly access the symmetry matrices (in lattice coords)
	const std::vector<matrix>& getSphericalMatrices(int l, bool relativistic) const; //!< directly access the symmetry matrices (in Ylm or spin-angle basis at specified l, depending on relativistic)
	const std::vector<int>& getKpointInvertList() const; //!< direct access to inversion property of symmetry group (see kpointInvertList)
	const std::vector<std::vector<std::vector<int> > >& getAtomMap() const; //!< direct access to mapping of each atom under each symmetry matrix (index order species, atom, symmetry)
	void printKmap(FILE* fp) const; //!< print the k-point map (cached in kmap)
	
	static matrix getSpinorRotation(const matrix3<>& rot); //calculate spinor rotation from Cartesian rotation matrix
private:
	const Everything* e;
	
	std::vector<SpaceGroupOp> sym; //!< set of space group operations
	std::vector< std::vector<matrix> > symSpherical; //!< symmetry matrices for real-Ylm basis objects
	std::vector< std::vector<matrix> > symSpinAngle; //!< symmetry matrices for spin-angle-function basis objects (relativistic version of symSpherical)
	
	std::vector<int> kpointInvertList; //!< Contains +1 for empty or inversion-containing symmetry group, and contains +1 and -1 otherwise
	std::vector<unsigned long long> kmap; //!< mapping from unreduced mesh to reduced set under symmetries
	friend struct CommandSymmetries;
	friend struct CommandSymmetryMatrix;
	friend struct CommandDebug;
	friend struct CommandKpointReduceInversion;
	
	bool kReduceUseInversion; //!< whether to use inversion symmetry to reduce k-point mesh
	bool shouldPrintMatrices;
	
	void calcSymmetries(); //!< Calculate symmetries of the entire system
	
	//! Find space group for the lattice with basis, given the point group symmetries of lattice alone
	std::vector<SpaceGroupOp> findSpaceGroup(const std::vector< matrix3<int> >& symLattice) const; 

	void sortSymmetries(); //!< ensure that the first symmetry is identity
	void checkFFTbox(); //!< verify that the sampled mesh is commensurate with symmetries
	void checkSymmetries(); //!< check validity of manually specified symmetry matrices
	
	//Index map for scalar field (electron density, potential) symmetrization in reciprocal space
	IndexArray symmIndex; //sets of nSym consecutive G-space indices that should be averaged during symmetrization
	IndexArray symmMult; //multiplicity (how many times each element is repeated) in each equivalence class
	ManagedArray<complex> symmIndexPhase; //phase factor for entry at each index
	ManagedArray<matrix3<>> symmRotSpin; //nSym Cartesian (pseudo-vector) rotation matrices for spin-density symmetrization
	void initSymmIndex();
	
	//Atom maps:
	std::vector<std::vector<std::vector<int> > > atomMap;
	void initAtomMaps();
	
	//Supercell handling (for phonon):
	vector3<int> sup; //!< this is an exact supercell of some unit cell with this count and restrict space group to translations within that unit cell
	bool isPertSup; //!< whether this is a perturbed supercell; if so, manual symmetries are reduced by atom positions (instead of raising an error)
	friend class Phonon;
};

//! @}
#endif // JDFTX_ELECTRONIC_SYMMETRIES_H
