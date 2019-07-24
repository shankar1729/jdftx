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
struct WannierGradient
{	std::vector<matrix> B1;
	std::vector<matrix> B2;
	const class WannierMinimizer* wmin;
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
		matrix U, Omega_UdotU; //!< net rotation (nBands x nIn) and intermediate gradient (Omega_U * U) w.r.t it
		matrix U1; //Initial rotation into Wannier subspace (nBands x nIn)
		matrix U2; //Rotation within Wannier subspace (nIn x nIn)
		//NOTE: U and U2 are truncated from nIn -> nCenters after minimization
		std::shared_ptr<MPIUtil> mpi; //MPI communicator with head=whose(ik) over which U must be bcast'd and reduced (created by subclass if needed)
	};
	void bcastU(); //broadcast U to any processes that need it
	
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
	
	//Phonon related:
	std::vector<vector3<>> xAtoms;  //lattice coordinates of all atoms in order
	std::map<vector3<int>,matrix> ePhCellMap; //cell map for e-ph matrix elements
	int nPhononModes; //!< number of phonon modes
	diagMatrix invsqrtM; //!< 1/sqrt(M) per nuclear displacement mode
	int prodPhononSup; //number of unit cells in phonon supercell
	struct KpointPair { vector3<> k1, k2; int ik1, ik2; };
	std::vector<KpointPair> kpointPairs; //pairs of k-points in the same order as matrices in phononHsub
	int iPairStart, iPairStop; //MPI division for working on kpoint pairs during phonon processing
	
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
	
	//! Overlap between columnbundles of different k-points, with appropriate ultrasoft augmentation.
	//! Note that the augmentation in the O() from electronic/ColumnBundle.h assumes both sides have same k-point.
	//! If provided, use the cached projections instead of recomputing them.
	matrix overlap(const ColumnBundle& C1, const ColumnBundle& C2, const std::vector<matrix>* VdagC1ptr=0, const std::vector<matrix>* VdagC2ptr=0) const;
	
	//! Preconditioner for Wannier optimization: identity by default, override in derived class to change
	virtual WannierGradient precondition(const WannierGradient& grad);

private:

	//! Get the trial wavefunctions (hydrogenic, atomic or numerical orbitals) for the group of centers in the common basis
	ColumnBundle trialWfns(const Kpoint& kpoint) const;
	std::map< Kpoint, std::shared_ptr<ColumnBundle> > numericalOrbitals; //!< numerical orbitals read from file
	
	//! Return an exactly unitary version of U (orthogonalize columns)
	//! If isSingular is provided, function will set it to true and return rather than stack-tracing in singular cases.
	static matrix fixUnitary(const matrix& U, bool* isSingular=0); 
	
	//! Load / compute rotations for a given spin channel (used by saveMLWF)
	void initRotations(int iSpin);
	
	//! Wannierize and dump a Bloch-space matrix to file, optionally zeroing out the real parts
	void dumpWannierized(const matrix& Htilde, const matrix& phase, string varName, bool realPartOnly, int iSpin) const;
	
	//---- Shared variables and subroutines implementing various Wannier outputs within saveMLWF() ----
	bool realPartOnly; //whether outputs should have only real part
	std::vector<vector3<>> xExpect; //converged wannier centers in lattice coordinates
	std::map<vector3<int>,matrix> iCellMap; //cell map for wannier output accounting for xExpect
	std::vector<matrix> DblochMesh; //gradient matix elements in Bloch basis
	void saveMLWF_C(int iSpin); //Wavefunctions
	void saveMLWF_H(int iSpin, const matrix& phase); //Hamiltonian
	void saveMLWF_P(int iSpin, const matrix& phase); //Momenta
	void saveMLWF_D(int iSpin, const matrix& phase); //Gradient
	void saveMLWF_S(int iSpin, const matrix& phase); //Spins
	void saveMLWF_W(int iSpin, const matrix& phase); //Slab weights
	void saveMLWF_ImSigma_ee(int iSpin, const matrix& phase); //e-e linewidths
	void saveMLWF_phonon(int iSpin); //e-ph matrix elements and related
};


//! Long range sum over G-vectors used for polar corrections
class LongRangeSum
{
	const matrix3<> G, GGT, epsInfLat;
	const double alphaInvBy4;
	vector3<int> iGbox;
public:
	LongRangeSum(const matrix3<>& R, const matrix3<>& epsInf, const double alpha=0.1)
	: G(2*M_PI*inv(R)), GGT(G*(~G)),
		epsInfLat(G * epsInf * (~G)),
		alphaInvBy4(0.25/alpha)
	{	//Initialize sample counts
		const double Gmax = 12.*sqrt(alpha); //such that exp(-Gmax^2/(4 alpha)) < 1e-16
		for(int i=0; i<3; i++)
			iGbox[i] = int(Gmax*R.column(i).length()/(2*M_PI)) + 2; //margin to account for q
	}
	
	inline complex operator()(const vector3<>& q, const vector3<>& Zeff, const vector3<>& atpos)
	{	//Wrap q to fundamental zone
		vector3<> qBZ = q;
		for(int iDir=0; iDir<3; iDir++)
			qBZ[iDir] -= floor(qBZ[iDir] + 0.5);
		//Convert Z to lattice coordinates:
		vector3<> ZeffLat = G * Zeff;
		//Loop over G-vectors:
		complex result = 0.;
		vector3<int> iG;
		for(iG[0]=-iGbox[0]; iG[0]<=iGbox[0]; iG[0]++)
		for(iG[1]=-iGbox[1]; iG[1]<=iGbox[1]; iG[1]++)
		for(iG[2]=-iGbox[2]; iG[2]<=iGbox[2]; iG[2]++)
		{	vector3<> iGq = iG + qBZ;
			double expFac = alphaInvBy4 * GGT.metric_length_squared(iGq);
			if(expFac > 1e-16 and expFac < 36.) //i.e. G != 0 and gaussian non-zero at double precision
				result += cis((-2*M_PI)*dot(iGq,atpos)) * exp(-expFac) * dot(ZeffLat, iGq) / epsInfLat.metric_length_squared(iGq);
		}
		return result;
	}
};

//! @}
#endif // JDFTX_WANNIER_WANNIERMINIMIZER_H
