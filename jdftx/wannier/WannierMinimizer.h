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

#ifndef JDFTX_WANNIER_WANNIERMINIMIZER_H
#define JDFTX_WANNIER_WANNIERMINIMIZER_H

#include <core/Util.h>
#include <core/Minimize.h>
#include <core/LatticeUtils.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/matrix.h>
#include <electronic/ColumnBundleTransform.h>
#include <wannier/Wannier.h>

//! @addtogroup Output
//! @{
//! @file WannierMinimizer.h

//! State / vector space entry for wannier minimization
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

//! Base class for different wannier minimizers:
class WannierMinimizer : public Minimizable<WannierGradient>
{
public:
	WannierMinimizer(const Everything& e, const Wannier& wannier, bool needSuperOverride=false);
	virtual ~WannierMinimizer() {}
	void initTransformDependent(); //!< second half of common initialization that must happen after sub-class is fully initialized
	
	void saveMLWF(); //!< save wannier functions for all spins
	void saveMLWF(int iSpin); //!< save for specified spin
	
	//Interface for minimize:
	void step(const WannierGradient& grad, double alpha);
	double compute(WannierGradient* grad, WannierGradient* Kgrad); //!< identity preconditioner, but impose hermiticity constraint
	void constrain(WannierGradient& grad); //!< enforce hermiticity
	bool report(int iter);
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	//! Entries in the k-point mesh
	struct Kpoint : public QuantumNumber, public Supercell::KmeshTransform
	{	bool operator<(const Kpoint& other) const;
		bool operator==(const Kpoint& other) const;
	};
	
	//! Entry in the k-point mesh, including state of minimizer (subspace rotations)
	struct KmeshEntry
	{	Kpoint point; //!< point in the mesh
		//State of system for Wannier minimize:
		int nIn; //number of bands that contribute to the Wannier subspace
		int nFixed; //number of bands that contribute fully to the Wannier subspace (cannot be partially mixed out)
		matrix U, Omega_U; //!< net rotation (nBands x nCenters) and intermediate gradient w.r.t it
		matrix U1; //Initial rotation into Wannier subspace (nBands x nIn)
		matrix U2; //Rotation within Wannier subspace (nCenters x nCenters)
	};

	//-------- Interface for subclasses that provide the objective function for Wannier minimization
	
	//! Prepare for minimization of spin channel iSpin
	virtual void initialize(int iSpin)=0;
	
	//! Return Omega and set rExpect and rSqExpect.
	//! If grad=true, set Omega_U in the KMeshEntry array.
	//! Note base class handles computing U and propagating gradients.
	//! At input, U will be available on all processes and Omega_U will be zero.
	//! Base class will accumulate Omega_U across processes on return.
	virtual double getOmega(bool grad=false)=0;
	
	//! Like getOmega, but for the subspace-invariant OmegaI instead.
	virtual double getOmegaI(bool grad=false)=0;
	
protected:
	friend struct WannierGradient;
	const Everything& e;
	const Wannier& wannier;
	const std::vector<SpaceGroupOp>& sym;
	int nCenters, nFrozen, nBands; //!< number of Wannier centers (total and frozen) and source bands
	int nSpins, qCount; //!< number of spins, and number of states per spin
	int nSpinor; //!< number of spinor components
	std::vector<double> rSqExpect; //!< Expectation values for r^2 per center in current group
	std::vector< vector3<> > rExpect; //!< Expectation values for r per center in current group
	std::vector<bool> pinned; //!< Whether centers are pinned or free
	std::vector< vector3<> > rPinned; //!< Where centers are pinned to (if they are)
	
	//Supercell grid and basis:
	bool needSuper; //!< whether supercell is necessary
	GridInfo gInfoSuper; //!< supercell grid
	Basis basisSuper; //!< supercell wavefcuntion basis
	QuantumNumber qnumSuper; //!< supercell k-point
	int nPhononModes; //!< number of phonon modes
	diagMatrix invsqrtM; //!< 1/sqrt(M) per nuclear displacement mode
	
	std::vector<KmeshEntry> kMesh; //!< k-point mesh with FD formula
	std::set<Kpoint> kpoints; //!< list of all k-points that will be in use (including those in FD formulae)
	std::shared_ptr<ColumnBundleTransform::BasisWrapper> basisWrapper, basisSuperWrapper; //!< look-up tables for initializing transforms
	std::map<Kpoint, std::shared_ptr<ColumnBundleTransform> > transformMap, transformMapSuper; //!< wave-function transforms for each k-point to the common bases
	Basis basis; //!< common basis (with indexing into full G-space)
	
	//k-mesh MPI division:
	TaskDivision kDivision;
	size_t ikStart, ikStop;
	inline bool isMine(size_t ik) const { return kDivision.isMine(ik); }
	inline int whose(size_t ik) const { return kDivision.whose(ik); }

	//state MPI division (wrappers to ElecInfo)
	bool isMine_q(int ik, int iSpin) const { return e.eInfo.isMine(kMesh[ik].point.iReduced + iSpin*qCount); }
	int whose_q(int ik, int iSpin) const { return e.eInfo.whose(kMesh[ik].point.iReduced + iSpin*qCount); }
	
	//! Get the wavefunctions for a particular k-point in the common basis (and optionally retrieve psp projections)
	ColumnBundle getWfns(const Kpoint& kpoint, int iSpin, std::vector<matrix>* VdagResultPtr=0) const;
	std::vector<ColumnBundle> Cother; //!< wavefunctions from another process
	std::vector<std::vector<matrix>> VdagCother; //!< psp projections of wavefunctions from another process
	
	//! Like getWfns, but accumulate instead of setting, and with optional transformation matrix: result += alpha * wfns * A
	void axpyWfns(double alpha, const matrix& A, const Kpoint& kpoint, int iSpin, ColumnBundle& result, std::vector<matrix>* VdagResultPtr=0) const;
	
	//! Gradient propagation corresponding to axpyWfns: from dOmega/d(result) to dOmega/dA
	void axpyWfns_grad(double alpha, matrix& Omega_A, const Kpoint& kpoint, int iSpin, const ColumnBundle& Omega_result) const;
	
	//! Get the trial wavefunctions (hydrogenic, atomic or numerical orbitals) for the group of centers in the common basis
	ColumnBundle trialWfns(const Kpoint& kpoint) const;
	std::map< Kpoint, std::shared_ptr<ColumnBundle> > numericalOrbitals; //!< numerical orbitals read from file
	
	//! Overlap between columnbundles of different k-points, with appropriate ultrasoft augmentation.
	//! Note that the augmentation in the O() from electronic/ColumnBundle.h assumes both sides have same k-point.
	//! If provided, use the cached projections instead of recomputing them.
	matrix overlap(const ColumnBundle& C1, const ColumnBundle& C2, const std::vector<matrix>* VdagC1ptr=0, const std::vector<matrix>* VdagC2ptr=0) const;
	
	//! Dump a named matrix variable to file, optionally zeroing out the real parts
	void dumpMatrix(const matrix& H, string varName, bool realPartOnly, int iSpin) const;
	
	//! Return an exactly unitary version of U (orthogonalize columns)
	//! If isSingular is provided, function will set it to true and return rather than stack-tracing in singular cases.
	static matrix fixUnitary(const matrix& U, bool* isSingular=0); 
	
	//! Preconditioner for Wannier optimization: identity by default, override in derived class to change
	virtual WannierGradient precondition(const WannierGradient& grad);
};

//! @}
#endif // JDFTX_WANNIER_WANNIERMINIMIZER_H
