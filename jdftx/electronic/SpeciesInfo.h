/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler
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

#ifndef JDFTX_ELECTRONIC_SPECIESINFO_H
#define JDFTX_ELECTRONIC_SPECIESINFO_H

#include <electronic/VanDerWaals.h>
#include <core/RadialFunction.h>
#include <core/matrix.h>
#include <core/ScalarFieldArray.h>
#include <core/vector3.h>
#include <core/string.h>

class ColumnBundle;
class QuantumNumber;
class Basis;

//! @addtogroup IonicSystem
//! @{

//! Pseudopotential for a species of ions, and atom positions and other properties for that species
class SpeciesInfo
{
public:

	double Z; //!< Valence charge of the species (prefactor to 1/r in long-range part of pseudopotential)
	int atomicNumber; //!< Atomic number of the species (0 if unavailable)
	string name; //!< Identifier
	string potfilename; //!< pseudopotential filename
	string pulayfilename; //!< pulay correction filename
	bool fromWildcard; //!< whether this pseudopotential was automatically added using a wildcard (for command printing purposes only)
	
	std::vector<vector3<> > atpos; //!< array of atomic positions of this species
	std::vector<vector3<> > velocities; //!< array of atomic velocities (null unless running MD) in lattice coordinates
	ManagedArray<vector3<>> atposManaged; //!< managed copy of atpos accessed from operator code (for auto cpu/gpu transfers)
	void sync_atpos(); //!< update changes in atpos; call whenever atpos is changed (this will update atposManaged and invalidate cached projectors, if any)
	
	double dE_dnG; //!< Derivative of [total energy per atom] w.r.t [nPlanewaves per unit volume] (for Pulay corrections)
	double mass; //!< ionic mass (currently unused)	
	double coreRadius; //!< maximum pseudopotential core radius (used for core overlap checks during ionic/lattice relaxation)
	double ZfullCore; //!< number of electrons in full-core correction (atomicNumber - Z - integral(nCore))
	
	//! Contains the information on the constraints of motion for each ion.
	struct Constraint
	{	double moveScale; //! preconditioning factor (0 fixes ion)
		vector3<> d; //! The direction or plane normal of the constraint (in cartesian coordinates)
		enum ConstraintType {None, Linear, Planar, HyperPlane} type; //! Type of the constraint that is being applied to the ion
		
		vector3<> operator()(const vector3<>& grad) const;
		bool isEquivalent(const Constraint& otherConstraint, const matrix3<>& transform) const;
		int getDimension() const;
		void print(FILE* fp, const Everything&) const;
	};
	std::vector<Constraint> constraints; //!< List of all constraints on ions of this species
	
	std::vector< vector3<> > initialMagneticMoments; //!< Initial magnetic moments of each atom (used only for LCAO and symmetries) (x and y magnetizations only used in noncollinear calculations)
	double initialOxidationState; //!< Initial oxidation state of this species (only affects LCAO)
	
	SpeciesInfo();
	~SpeciesInfo();
	void setup(const Everything&);
	void print(FILE* fp) const; //!< print ionic positions from current species
	void populationAnalysis(const std::vector<matrix>& RhoAll) const; //!< print population analysis given the density matrix in the Lowdin basis
	bool isRelativistic() const { return psi2j.size(); } //!< whether pseudopotential is relativistic
	bool isUltrasoft() const { return Qint.size(); } //!< whether pseudopotential is ultrasoft
	
	enum PseudopotentialFormat
	{	Fhi, //!< FHI format with ABINIT header (.fhi files)
		Uspp, //!< USPP format ultrasoft pseudopotential
		UPF //!< Quantum Espresso's Universal Pseudopotential Format (the XML version 2 only)
	};
	//! Returns the pseudopotential format
	PseudopotentialFormat getPSPFormat(){return pspFormat;}

	std::shared_ptr<ColumnBundle> getV(const ColumnBundle& Cq) const; //!< get projectors with qnum and basis matching Cq  (optionally cached)
	int nProjectors() const { return MnlAll.nRows() * atpos.size(); } //!< total number of projectors for all atoms in this species (number of columns in result of getV)
	
	//! Return non-local energy for this species and quantum number q and optionally accumulate
	//! projected electronic gradient in HVdagCq (if non-null)
	double EnlAndGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq, matrix& HVdagCq) const;
	
	//! Accumulate pseudopotential contribution to the overlap in OCq
	void augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq, matrix* VdagCq=0) const;
	
	//! Clear internal data and prepare for density augmentation (call before a loop over augmentDensitySpherical per k-point)
	void augmentDensityInit();
	//! Accumulate the pseudopotential dependent contribution to the density in the spherical functions nAug (call once per k-point)
	void augmentDensitySpherical(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq);
	//! Accumulate the spherical augmentation functions nAug to the grid electron density (call only once, after augmentDensitySpherical on all k-points)
	void augmentDensityGrid(ScalarFieldArray& n) const;
	
	//! Gradient propagation corresponding to augmentDensityGrid (stores intermediate spherical function results to E_nAug; call only once) 
	void augmentDensityGridGrad(const ScalarFieldArray& E_n, std::vector<vector3<> >* forces=0);
	//! Gradient propagation corresponding to augmentDensitySpherical (uses intermediate spherical function results from E_nAug; call once per k-point after augmentDensityGridGrad) 
	void augmentDensitySphericalGrad(const QuantumNumber& qnum, const matrix& VdagCq, matrix& HVdagCq) const;
	
	//DFT+U functions: handle IonInfo::rhoAtom_*() for this species
	//The rhoAtom pointers point to the start of those relevant to this species (and ends at that pointer + rhoAtom_nMatrices())
	size_t rhoAtom_nMatrices() const;
	void rhoAtom_initZero(matrix* rhoAtomPtr) const;
	void rhoAtom_calc(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, matrix* rhoAtomPtr) const;
	double rhoAtom_computeU(const matrix* rhoAtomPtr, matrix* U_rhoAtomPtr) const;
	void rhoAtom_grad(const ColumnBundle& Cq, const matrix* U_rhoAtomPtr, ColumnBundle& HCq) const;
	void rhoAtom_forces(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, const matrix* U_rhoAtomPtr, std::vector<vector3<> >& forces) const;
	void rhoAtom_getV(const ColumnBundle& Cq, const matrix* U_rhoAtomPtr, ColumnBundle& Opsi, matrix& M) const; //get DFT+U Hamiltonian in the same format as the nonlocal pseudopotential (psi = atomic orbitals, M = matrix in that order)

	//Atomic orbital related functions:
	void accumulateAtomicDensity(ScalarFieldTildeArray& nTilde) const; //!< Accumulate atomic density from this species
	void accumulateAtomicPotential(ScalarFieldTilde& dTilde) const; //!< Accumulate electrostatic potential of neutral atoms from this species
	void setAtomicOrbitals(ColumnBundle& Y, bool applyO, int colOffset=0) const; //!< Calculate atomic orbitals with/without O preapplied (store in Y with an optional column offset)
	void setAtomicOrbitals(ColumnBundle& Y, bool applyO, unsigned n, int l, int colOffset=0, int atomColStride=0) const;  //!< Same as above, but for specific n and l.
		//!< If non-zero, atomColStride overrides the number of columns between the same orbital of multiple atoms (default = number of orbitals at current n and l)
	int nAtomicOrbitals() const; //!< return number of atomic orbitals in this species (all atoms)
	int lMaxAtomicOrbitals() const; //!< return maximum angular momentum in available atomic orbitals
	int nAtomicOrbitals(int l) const; //!< return number of (pseudo-)principal quantum numbers for atomic orbitals of given l
	int atomicOrbitalOffset(unsigned iAtom, unsigned n, int l, int m, int s) const; //!< offset of specified atomic orbital in output of current species (when not using the fixed n and l version)
		//!< s is 0/1 for up/dn spinors in non-relativistic case, s=0/1 is for j=l+/-0.5 and mj=m+/-0.5 in relativistic case

	//! Add contributions from this species to Vlocps, rhoIon, nChargeball and nCore/tauCore (if any)
	void updateLocal(ScalarFieldTilde& Vlocps, ScalarFieldTilde& rhoIon, ScalarFieldTilde& nChargeball,
		ScalarFieldTilde& nCore, ScalarFieldTilde& tauCore) const; 
	
	//! Return the local forces (due to Vlocps, rhoIon, nChargeball and nCore/tauCore)
	std::vector< vector3<> > getLocalForces(const ScalarFieldTilde& ccgrad_Vlocps, const ScalarFieldTilde& ccgrad_rhoIon,
		const ScalarFieldTilde& ccgrad_nChargeball, const ScalarFieldTilde& ccgrad_nCore, const ScalarFieldTilde& ccgrad_tauCore) const;

	//! Propagate gradient with respect to atomic projections (in E_VdagC, along with additional overlap contributions from grad_CdagOC) to forces:
	void accumNonlocalForces(const ColumnBundle& Cq, const matrix& VdagC, const matrix& E_VdagC, const matrix& grad_CdagOCq, std::vector<vector3<> >& forces) const;
	
	//! Spin-angle helper functions:
	static matrix getYlmToSpinAngleMatrix(int l, int j2); //!< Get the ((2l+1)*2)x(j2+1) matrix that transforms the Ylm+spin to the spin-angle functions, where j2=2*j with j = l+/-0.5
	static matrix getYlmOverlapMatrix(int l, int j2); //!< Get the ((2l+1)*2)x((2l+1)*2) overlap matrix of the spin-spherical harmonics for total angular momentum j (note j2=2*j)
private:
	matrix3<> Rprev; void updateLatticeDependent(); //!< If Rprev differs from gInfo.R, update the lattice dependent quantities (such as the radial functions)

	RadialFunctionG VlocRadial; //!< local pseudopotential
	RadialFunctionG nCoreRadial; //!< core density for partial core correction
	RadialFunctionG tauCoreRadial; //!< core KE density for partial core correction with meta-GGAs

	std::vector< std::vector<RadialFunctionG> > VnlRadial; //!< non-local projectors (outer index l, inner index projetcor)
	std::vector<matrix> Mnl; //!< nonlocal pseudopotential projector matrix (indexed by l)
	matrix MnlAll; //!< block matrix containing Mnl for all l,m 
	
	std::vector<matrix> Qint; //!< overlap augmentation matrix (indexed by l, empty if no augmentation)
	matrix QintAll; //!< block matrix containing Qint for all l,m 
	
	std::map<std::pair<vector3<>,const Basis*>, std::shared_ptr<ColumnBundle> > cachedV; //cached projectors (identified by k-point and basis pointer)
	
	struct QijIndex
	{	int l1, p1; //!< Angular momentum and projector index for channel i
		int l2, p2; //!< Angular momentum and projector index for channel j
		int l; //!< net angular momentum (l varies from |l1-l2| to (l1+l2) in increments of 2)
		int index; //!< psoition in Qradial map
		bool operator<(const QijIndex&) const; //!< comparison that ensures (l1,p1)<-->(l2,p2) symmetry
	private:
		void sortIndices(); //!< swap (l1,p1)<-->(l2,p2) indices if needed to bring to the upper triangular part
	};
	std::map<QijIndex,RadialFunctionG> Qradial; //!< radial functions for density augmentation
	matrix QradialMat; //!< matrix with all the radial augmentation functions in columns (ordered by index)
	matrix nAug; //!< intermediate electron density augmentation in the basis of Qradial functions (Flat array indexed by spin, atom number and then Qradial index)
	matrix E_nAug; //!< Gradient w.r.t nAug (same layout)
	ManagedArray<uint64_t> nagIndex; ManagedArray<size_t> nagIndexPtr; //!< grid indices arranged by |G|, used for coordinating scattered accumulate in nAugmentGrad(_gpu)

	std::vector<std::vector<RadialFunctionG> > psiRadial; //!< radial part of the atomic orbitals (outer index l, inner index shell)
	std::vector<std::vector<RadialFunctionG> > OpsiRadial; //!< O(psiRadial): includes Q contributions for ultrasoft pseudopotentials
	std::vector<std::vector<double> > atomEigs; //!< Eigenvalues of the atomic orbitals in the atomic state (read in, or computed from tail of psi by estimateAtomEigs)
	
	//! Extra information for spin-orbit coupling:
	std::vector< std::vector<int> > Vnl2j, psi2j; //!< 2*j for the projectors and orbitals respectively
	matrix fljAll; //!< block-diagonal collection of YlmOverlapMatrix for all projectors on one atom
	
	//!Parameters for optional DFT+U corrections
	struct PlusU
	{	int n, l; //!< Principal and angular quantum numbers
		double UminusJ; //!< U - J [Dudarev et al, Phys. Rev. B 57, 1505] scheme
		
		std::vector<double> Vext; //! external potential per atom (used for calculating linear-response U)
	};
	std::vector<PlusU> plusU; //!< list of +U corrections
	
	std::shared_ptr<VanDerWaals::AtomParams> vdwOverride; //!< optional override of vdW parameters for this species
	std::shared_ptr<double> atomicRadiusOverride; //!< optional override of atomic radius (used for some solvation models) for this species
	
	PseudopotentialFormat pspFormat;
	
	// gaussian chargeball used to prevent dielectric spikes
	// DEPRECATED: As far as possible, use partial core correction instead
	double Z_chargeball; //!< normalization of chargeball
	double width_chargeball; //!< width of gaussian chargeball

	const Everything* e;
	
	double tauCore_rCut; bool tauCorePlot; //!< radial cutoff for generating core kinetic energy density and whether to generate a plot file of the result
	
	void setCore(RadialFunctionR&); //!< Generate nCoreRadial and tauCoreRadial if required (in SpeciesInfo_core.cpp)
	void readFhi(istream&); // Implemented in SpeciesInfo_readFhi.cpp
	void readUspp(istream&); //Implemented in SpeciesInfo_readUspp.cpp
	void readUPF(istream&); //Implemented in SpeciesInfo_readUPF.cpp
	void setupPulay();
	void setPsi(std::vector<std::vector<RadialFunctionR> >& psi); //!< Normalize, transform from real space and set psiRadial, OpsiRadial (for ultrasoft psps, call after setting Qint)
	
	//Following implemented in SpeciesInfo_atomFillings.cpp
	void estimateAtomEigs(); //!< If not read from file, estimate atomic eigenvalues from orbitals.
	void getAtom_nRadial(int spin, double magneticMoment, RadialFunctionG& nRadial, bool forceNeutral) const; //!< Compute the atomic density per spin channel, given the magnetic moment
	void getAtomPotential(RadialFunctionG& dRadial) const; //!< Get the total electrostatic potential of a neutral atom

	friend struct CommandAddU;
	friend struct CommandChargeball;
	friend struct CommandIonSpecies;
	friend struct CommandSetVDW;
	friend struct CommandSetAtomicRadius;
	friend struct CommandTauCore;
	friend struct CommandWavefunction;
	friend struct ElectronScattering;
	friend struct FluidSolverParams;
	friend class ColumnBundleTransform;
	friend class Dump;
	friend class BGW;
	friend class IonicMinimizer;
	friend class IonInfo;
	friend class PCM;
	friend class Phonon;
	friend class VanDerWaals;
	friend class WannierMinimizer;
};


static EnumStringMap<SpeciesInfo::Constraint::ConstraintType> constraintTypeMap
(	SpeciesInfo::Constraint::None, "None",
	SpeciesInfo::Constraint::Linear, "Linear",
	SpeciesInfo::Constraint::Planar, "Planar",
	SpeciesInfo::Constraint::HyperPlane, "HyperPlane"
);

//! @}
#endif // JDFTX_ELECTRONIC_SPECIESINFO_H
