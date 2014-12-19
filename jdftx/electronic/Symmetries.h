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

#include <electronic/common.h>
#include <core/matrix3.h>
#include <core/Data.h>
#include <vector>
#include <list>

enum SymmetryMode {SymmetriesNone, SymmetriesAutomatic, SymmetriesManual}; //!< symmetry modes

class Symmetries
{
public:
	SymmetryMode mode; //!< Symmetry mode (none, automatic or manual)

	Symmetries();
	~Symmetries();
	void setup(const Everything& everything); //!< Phase 1 of setup which computes/checks lattice+basis symmetries
	void setupMesh(); //!< Phase 2 of setup which computes / checks FFTbox and k-mesh dependent symmetries
	
	//! Reduce a k-point mesh (and remember its inversion symmetry property in kpointInvertList)
	std::vector<QuantumNumber> reduceKmesh(const std::vector<QuantumNumber>& qnums) const;
	
	void symmetrize(DataRptr&) const; //!< symmetrize a scalar field
	void symmetrize(IonicGradient&) const; //!< symmetrize forces
	void symmetrizeSpherical(matrix&, const SpeciesInfo* specie) const; //!< symmetrize matrices in Ylm basis per atom of species sp (accounting for atom maps)
	const std::vector< matrix3<int> >& getMatrices() const; //!< directly access the symmetry matrices (in lattice coords)
	const std::vector< matrix3<int> >& getMeshMatrices() const; //!< directly access the symmetry matrices (in mesh coords)
	const std::vector<matrix>& getSphericalMatrices(int l, bool relativistic) const; //!< directly access the symmetry matrices (in Ylm or spin-angle basis at specified l, depending on relativistic)
	const std::vector<int>& getKpointInvertList() const; //!< direct access to inversion property of symmetry group (see kpointInvertList)
	const std::vector<std::vector<std::vector<int> > >& getAtomMap() const; //!< direct access to mapping of each atom under each symmetry matrix (index order species, atom, symmetry)
	void printKmap(FILE* fp) const; //!< print the k-point map (cached in kmap)
	
	static matrix getSpinorRotation(const matrix3<>& rot); //calculate spinor rotation from Cartesian rotation matrix
private:
	const Everything* e;
	std::vector< matrix3<int> > sym; //!< symmetry matrices in covariant lattice coordinates
	std::vector< matrix3<int> > symMesh; //!< symmetry matrices in covariant mesh coordinates (lattice / sample counts)
	std::vector< std::vector<matrix> > symSpherical; //!< symmetry matrices for real-Ylm basis objects
	std::vector< std::vector<matrix> > symSpinAngle; //!< symmetry matrices for spin-angle-function basis objects (relativistic version of symSpherical)
	
	std::vector<int> kpointInvertList; //!< Contains +1 for empty or inversion-containing symmetry group, and contains +1 and -1 otherwise
	std::vector<unsigned long long> kmap; //!< mapping from unreduced mesh to reduced set under symmetries
	friend class CommandSymmetries;
	friend class CommandSymmetryMatrix;
	friend class CommandDebug;
	
	bool shouldPrintMatrices;
	bool shouldMoveAtoms;
	
	void calcSymmetries(); //!< Calculate symmetries of the entire system
	
	//! Find subgroup of lattice symmetries for the lattice with basis (optionally offset by some amount)
	std::vector< matrix3<int> > basisReduce(const std::vector< matrix3<int> >& symLattice, vector3<> offset=vector3<>()) const; 

	void sortSymmetries(); //!< ensure that the first symmetry is identity
	void checkFFTbox(); //!< verify that the sampled mesh is commensurate with symmetries (and generate symMesh)
	void checkSymmetries() const; //!< check validity of manually specified symmetry matrices
	
	//Index map for scalar field (electron density, potential) symmetrization
	int *symmIndex, nSymmIndex;
	void initSymmIndex();
	
	//Atom maps:
	std::vector<std::vector<std::vector<int> > > atomMap;
	void initAtomMaps();
};

#endif // JDFTX_ELECTRONIC_SYMMETRIES_H
