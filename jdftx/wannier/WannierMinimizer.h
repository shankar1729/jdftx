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

#ifndef JDFTX_ELECTRONIC_WANNIERMINIMIZER_H
#define JDFTX_ELECTRONIC_WANNIERMINIMIZER_H

#include <core/Util.h>
#include <core/Minimize.h>
#include <core/LatticeUtils.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <wannier/Wannier.h>

struct WannierGradient : public std::vector<matrix>
{	const class WannierMinimizer* wmin;
    void init(const class WannierMinimizer* wmin);
	size_t ikStart() const;
	size_t ikStop() const;
};

//---- linear algebra functions required by Minimizable<WannierGradient> -----
WannierGradient clone(const WannierGradient& grad);
double dot(const WannierGradient& x, const WannierGradient& y);
WannierGradient& operator*=(WannierGradient& x, double alpha);
void axpy(double alpha, const WannierGradient& x, WannierGradient& y);
matrix randomMatrix(int nRows, int nCols);
void randomize(WannierGradient& x);

class WannierMinimizer : public Minimizable<WannierGradient>
{
public:
	WannierMinimizer(const Everything& e, const Wannier& wannier);
	virtual ~WannierMinimizer() {}
	
	void saveMLWF(); //save wannier functions for all spins
	void saveMLWF(int iSpin); //save for specified spin
	
	//Interface for minimize:
	void step(const WannierGradient& grad, double alpha);
	double compute(WannierGradient* grad);
	double sync(double x) const { mpiUtil->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	//! Entries in the k-point mesh
	struct Kpoint : public QuantumNumber, public Supercell::KmeshTransform
	{	bool operator<(const Kpoint& other) const;
		bool operator==(const Kpoint& other) const;
	};
	
private:
	friend class WannierGradient;
	const Everything& e;
	const Wannier& wannier;
	const std::vector< matrix3<int> >& sym;
	int nCenters, nBands; //!< number of Wannier centers and source bands
	int nSpins, qCount; //!< number of spins, and number of states per spin
	std::vector<double> rSqExpect; //!< Expectation values for r^2 per center in current group
	std::vector< vector3<> > rExpect; //!< Expectation values for r per center in current group
	double OmegaI; //invariant part
	
	//Supercell grid and basis:
	bool needSuper;
	GridInfo gInfoSuper;
	Basis basisSuper;
	QuantumNumber qnumSuper;
	std::map<vector3<int>,double> iCellMap; //unit-cell indices in supercell (and weights to account for surface multiplicity)
	
	//!An edge of the k-mesh involved in the finite difference formula
	struct EdgeFD
	{	double wb; //!< weight of neighbour
		vector3<> b; //!< displacement to neighbour
		unsigned ik; //!< index of neighbour in kMesh
		Kpoint point; //!< description of neighbour (source state, rotation, translation etc.)
		matrix M0; //!< initial overlap matrix for this pair (after applying the initial unitary rotations)
	};
	
	//!Entries in the k-point mesh with FD formula
	struct KmeshEntry
	{	Kpoint point; //!< point in the mesh
		std::vector<EdgeFD> edge; //!< neighbours with weights defining finite difference formula
		//State of system for Wannier minimize:
		int nIn; //number of bands that contribute to the Wannier subspace
		int nFixed; //number of bands that contribute fully to the Wannier subspace (cannot be partially mixed out)
		matrix B; //!< Independent variable for minimization (nCenters x nIn)
		matrix U, Omega_U; //!< net rotation (nBands x nCenters) and intermediate gradient w.r.t it
		//Stage 1: Select linear cominations of bands that enter Wannier subspace
		matrix U1, V1, B1evecs; //U1 = initial rotation (nBands x nIn), V1 = subsequent rotation from B (nIn x nCenters)
		diagMatrix B1eigs;
		//Stage 2: Rotations within Wannier subspace (all nCenters x nCenters)
		matrix U2, V2, B2evecs;
		diagMatrix B2eigs;
	};
	
	//!Indices from reduced basis to full G-space or a union of all reduced bases
	//!The super-suffixed versions indicate indices into the G-sphere/fftbox of the supercell
	struct Index
	{	const int nIndices;
		int *data, *dataSuper;
		#ifdef GPU_ENABLED
		int *dataGpu, *dataSuperGpu;
		#endif
		int *dataPref, *dataSuperPref;
		
		Index(int nIndices=0, bool needSuper=false); //!< allocate space for indices (optionally for supercell version as well)
		~Index();
		void set(); //!< update gpu-copies
		
		//Non-copyable:
		Index(const Index&)=delete;
		Index& operator=(const Index&)=delete;
	};
	
	std::vector<KmeshEntry> kMesh; //!< k-point mesh with FD formula
	std::map<Kpoint, std::shared_ptr<Index> > indexMap; //!< wave-function indexing from each k-point to the common union basis
	Basis basis; //!< common basis (with indexing into full G-space)
	
	//k-mesh MPI division:
	size_t ikStart, ikStop;
	std::vector<size_t> ikStopArr;
	bool isMine(size_t ik) const { return ik>=ikStart && ik<ikStop; }
	int whose(size_t ik) const;

	//state MPI division (wrappers to ElecInfo)
	bool isMine_q(int ik, int iSpin) const { return e.eInfo.isMine(kMesh[ik].point.iReduced + iSpin*qCount); }
	int whose_q(int ik, int iSpin) const { return e.eInfo.whose(kMesh[ik].point.iReduced + iSpin*qCount); }
	
	void addIndex(const Kpoint& kpoint); //!< Add index for a given kpoint to indexMap, with indices pointing to full G-space
	
	//! Get the wavefunctions for a particular k-point in the common basis
	ColumnBundle getWfns(const Kpoint& kpoint, int iSpin) const;
	std::vector<ColumnBundle> Cother; //wavefunctions from another process
	
	 //Like getWfns, but accumulate instead of setting, and with optional transformation matrix: result += alpha * wfns * A
	void axpyWfns(double alpha, const matrix& A, const Kpoint& kpoint, int iSpin, ColumnBundle& result) const;
	
	//! Get the trial wavefunctions (hydrogenic, atomic or numerical orbitals) for the group of centers in the common basis
	ColumnBundle trialWfns(const Kpoint& kpoint) const;
	std::map< Kpoint, std::shared_ptr<ColumnBundle> > numericalOrbitals;
	
	//! Overlap between columnbundles of different k-points, with appropriate ultrasoft augmentation
	//! (Note that the augmentation in the O() from electronic/operators.h assumes both sides have same k-point)
	matrix overlap(const ColumnBundle& C1, const ColumnBundle& C2) const;
};

#endif // JDFTX_ELECTRONIC_WANNIERMINIMIZER_H
