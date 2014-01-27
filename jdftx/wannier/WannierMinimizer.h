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

typedef std::vector<matrix> WannierGradient;

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
	const Everything& e;
	const Wannier& wannier;
	const std::vector< matrix3<int> >& sym;
	int nCenters, nBands; //!< number of Wannier centers and source bands
	int nSpins, qCount; //!< number of spins, and number of states per spin
	std::vector<double> rSqExpect; //!< Expectation values for r^2 per center in current group
	std::vector< vector3<> > rExpect; //!< Expectation values for r per center in current group
	
	//Supercell grid and basis:
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
		matrix B; //!< Independent variable for minimization
		matrix Bevecs; //!< Eigenvectors of B
		diagMatrix Beigs; //!< Eigenvalues of B
		matrix U0; //!< Initial unitary rotation (from trial wave functions)
		matrix V; //!< Subsequent unitary rotation of k-point given by e^iB (net unitary rotation = U0 * V)
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
	bool isMine(size_t ik) const { return ik>=ikStart && ik<ikStop; }

	//state MPI division (wrappers to ElecInfo)
	bool isMine_q(int ik, int iSpin) const { return e.eInfo.isMine(kMesh[ik].point.iReduced + iSpin*qCount); }
	int whose_q(int ik, int iSpin) const { return e.eInfo.whose(kMesh[ik].point.iReduced + iSpin*qCount); }
	
	void addIndex(const Kpoint& kpoint); //!< Add index for a given kpoint to indexMap, with indices pointing to full G-space
	
	//! Get the wavefunctions for a particular k-point for bands involved in current group
	//! The wavefunctions are returned in the common basis by default and in the supercell basis if super=true
	ColumnBundle getWfns(const Kpoint& kpoint, int iSpin, bool super=false) const;
	std::vector<ColumnBundle> Cother; //wavefunctions from another process
	
	//! Get the trial wavefunctions (gaussians) for the group of centers in the common basis
	ColumnBundle trialWfns(const Kpoint& kpoint) const;
	
	//! Overlap between columnbundles of different k-points, with appropriate ultrasoft augmentation
	//! (Note that the augmentation in the O() from electronic/operators.h assumes both sides have same k-point)
	matrix overlap(const ColumnBundle& C1, const ColumnBundle& C2) const;
};

#endif // JDFTX_ELECTRONIC_WANNIERMINIMIZER_H
