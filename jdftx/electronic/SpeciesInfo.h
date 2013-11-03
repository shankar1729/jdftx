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

#include <electronic/common.h>
#include <electronic/RadialFunction.h>
#include <electronic/matrix.h>
#include <core/DataCollection.h>
#include <core/vector3.h>
#include <core/string.h>


class SpeciesInfo
{
public:

	double Z; //!< Valence charge of the species (prefactor to 1/r in long-range part of pseudopotential)
	int atomicNumber; //!< Atomic number of the species (0 if unavailable)
	string name; //!< Identifier
	string potfilename, pulayfilename;
	
	std::vector<vector3<> > atpos; //!< array of atomic positions of this species
	#ifdef GPU_ENABLED
	vector3<> *atposGpu; //!< copy of atomic positions on the gpu
	void sync_atposGpu(); //!< update changes in atpos to atposGpu (call whenever atpos is changed)
	#endif
	vector3<>* atposPref; //!< points to atposGpu in GPU mode and atpos otherwise
	
	double dE_dnG; //!< Derivative of [total energy per atom] w.r.t [nPlanewaves per unit volume] (for Pulay corrections)
	double mass; //!< ionic mass (currently unused)	
	double coreRadius; //!< maximum pseudopotential core radius (used for core overlap checks during ionic/lattice relaxation)
	
	//! Contains the information on the constraints of motion for each ion.
	struct Constraint
	{	double moveScale; //! preconditioning factor (0 fixes ion)
		vector3<> d; //! The direction or plane normal of the constraint (in cartesian coordinates)
		enum ConstraintType {None, Linear, Planar} type; //! Type of the constraint that is being applied to the ion
		
		vector3<> operator()(const vector3<>& grad);
		bool isEquivalent(const Constraint& otherConstraint, const matrix3<>& transform) const;
		int getDimension() const;
		void print(FILE* fp, const Everything&) const;
	};
	std::vector<Constraint> constraints; //!< List of all constraints on ions of this species
	
	std::vector<double> initialMagneticMoments; //!< Initial magnetic moments of each atom in number of electrons (used only for LCAO and symmetries)
	
	SpeciesInfo();
	~SpeciesInfo();
	void setup(const Everything&);
	void print(FILE* fp) const;
	
	enum PseudopotentialFormat
	{	Pot, //!< Old DFT++ format (.pot files)
		Cpi, //!< FHI98 format (.cpi files)
		Fhi, //!< FHI format with ABINIT header (.fhi files)
		Uspp //!< USPP format ultrasoft pseudopotential
	};
	//! Returns the pseudopotential format
	PseudopotentialFormat getPSPFormat(){return pspFormat;}

	//! Return non-local energy for this species and quantum number q an doptionally accumulate
	//! electronic gradient in HCq (if non-null) and ionic gradient in forces (if non-null)
	double EnlAndGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, std::vector<vector3<> >* forces) const;
	
	//! Accumulate pseudopotential contribution to the overlap in OCq
	void augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const;
	
	//! Clear internal data and prepare for density augmentation (call before a loop ober augmentDensitySpherical per k-point)
	void augmentDensityInit();
	void augmentDensityCleanup();
	//! Accumulate the pseudopotential dependent contribution to the density in the spherical functions nAug (call once per k-point)
	void augmentDensitySpherical(const diagMatrix& Fq, const ColumnBundle& Cq);
	//! Accumulate the spherical augmentation functions nAug to the grid electron density (call only once, after augmentDensitySpherical on all k-points)
	void augmentDensityGrid(DataRptrCollection& n) const;
	
	//! Gradient propagation corresponding to augmentDensityGrid (stores intermediate spherical function results to E_nAug; call only once) 
	void augmentDensityGridGrad(const DataRptrCollection& E_n, std::vector<vector3<> >* forces=0);
	//! Gradient propagation corresponding to augmentDensitySpherical (uses intermediate spherical function results from E_nAug; call once per k-point after augmentDensityGridGrad) 
	void augmentDensitySphericalGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq,
		std::vector<vector3<> >* forces=0, const matrix& gradCdagOCq=matrix()) const;
	
	//! Perform IonInfo::computeU() for this species
	double computeU(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
		std::vector<ColumnBundle>* HC = 0, std::vector<vector3<> >* forces=0) const;
	
	//!Pseudo-atom configuration (retrieved along with atomic orbitals using setAtomicOrbitals)
	struct AtomConfig
	{	unsigned n; //!< pseudo-atom principal quantum number
		int l; //!< angular momentum
		unsigned iAtom; //!< atom index (within all atoms of this species)
		double F; //!< occupation in the atomic state (summed over m)
	};

	//! Accumulate atomic density from this species
	void accumulateAtomicDensity(DataGptrCollection& nTilde) const;
	//! Calculate atomic orbitals (store in Y with an optional column offset, and optionally retrieve the pseudo-atom configuration)
	void setAtomicOrbitals(ColumnBundle& Y, int colOffset=0, std::vector<AtomConfig>* atomConfig=0) const;
	//! Store a single atomic orbital (iAtom'th atom, n'th shell of angular momentum l with specified m value) in col'th column of Y:
	void setAtomicOrbital(ColumnBundle& Y, int col, unsigned iAtom, unsigned n, int l, int m) const;
	int nAtomicOrbitals() const; //!< return number of atomic orbitals in this species (all atoms)
	int lMaxAtomicOrbitals() const; //!< return maximum angular momentum in available atomic orbitals
	int nAtomicOrbitals(int l) const; //!< return number of atomic orbitals of given l (per atom)
	int atomicOrbitalOffset(unsigned iAtom, unsigned n, int l, int m) const; //!< offset of specified atomic orbital in output of current species

	//! Add contributions from this species to Vlocps, rhoIon, nChargeball and nCore/tauCore (if any)
	void updateLocal(DataGptr& Vlocps, DataGptr& rhoIon, DataGptr& nChargeball,
		DataGptr& nCore, DataGptr& tauCore) const; 
	
	//! Return the local forces (due to Vlocps, rhoIon, nChargeball and nCore/tauCore)
	std::vector< vector3<> > getLocalForces(const DataGptr& ccgrad_Vlocps, const DataGptr& ccgrad_rhoIon,
		const DataGptr& ccgrad_nChargeball, const DataGptr& ccgrad_nCore, const DataGptr& ccgrad_tauCore) const;

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
	
	std::map< vector3<>, std::shared_ptr<ColumnBundle> > cachedV; //cached projectors (by k-point, assuming projectors are spin-independent)
	std::shared_ptr<ColumnBundle> getV(const ColumnBundle& Cq) const; //get projectors with qnum and basis matching Cq  (optionally cached)
	
	struct QijIndex
	{	int l1, p1; //!< Angular momentum and projector index for channel i
		int l2, p2; //!< Angular momentum and projector index for channel j
		int l; //!< net angular momentum (l varies from |l1-l2| to (l1+l2) in increments of 2)
		bool operator<(const QijIndex&) const; //!< comparison that ensures (l1,p1)<-->(l2,p2) symmetry
	private:
		void sortIndices(); //!< swap (l1,p1)<-->(l2,p2) indices if needed to bring to the upper triangular part
	};
	std::map<QijIndex,RadialFunctionG> Qradial; //!< radial functions for density augmentation
	double* nAug; //!< intermediate electron density augmentation in spherical functions (Flat array indexed by spin, atom number and then spline coeff)
	double* E_nAug; //!< Gradient w.r.t nAug (same layout)
	uint64_t* nagIndex; size_t* nagIndexPtr; //!< grid indices arranged by |G|, used for coordinating scattered accumulate in nAugmentGrad_gpu

	std::vector<std::vector<RadialFunctionG> > psiRadial; //!< radial part of the atomic orbitals (outer index l, inner index shell)
	std::vector<std::vector<RadialFunctionG> >* OpsiRadial; //!< O(psiRadial): includes Q contributions for ultrasoft pseudopotentials
	std::vector<std::vector<double> > atomEigs; //!< Eigenvalues of the atomic orbitals in the atomic state (read in, or computed from tail of psi by estimateAtomEigs)
	
	//! Calculate O(atomic orbitals) of specific n and l
	void setOpsi(ColumnBundle& Y, unsigned n, int l) const;
	
	//!Parameters for optional DFT+U corrections
	struct PlusU
	{	int n, l; //!< Principal and angular quantum numbers
		double UminusJ; //!< U - J [Dudarev et al, Phys. Rev. B 57, 1505] scheme
	};
	std::vector<PlusU> plusU; //!< list of +U corrections
	
	PseudopotentialFormat pspFormat;
	int lLocCpi; //!< local channel l for CPI files (Default: -1 => last channel in file)
	int recStartLen, recStopLen; //!< record marker lengths for fortran binary sequential file format (for uspp)
	double dGloc; //!< q-spacing for the local channel (default: 0.02)
	double dGnl; //!< q-spacing for the non-local projetcors (default: 0.02)
	
	// gaussian chargeball used to prevent dielectric spikes
	// DEPRECATED: As far as possible, use partial core correction instead
	double Z_chargeball; //!< normalization of chargeball
	double width_chargeball; //!< width of gaussian chargeball

	const Everything* e;
	
	void setCore(RadialFunctionR&); //!< Generate nCoreRadial and tauCoreRadial if required (in SpeciesInfo_core.cpp)
	void readPot(istream&); // Implemented in SpeciesInfo_readPot.cpp
	void readCpi(istream&); // Implemented in SpeciesInfo_readFhi.cpp
	void readFhi(istream&); // Implemented in SpeciesInfo_readFhi.cpp
	void readUspp(istream&); //Implemented in SpeciesInfo_readUspp.cpp
	void setupPulay(); // Implemented in SpeciesInfo_readPot.cpp
	
	//Following implemented in SpeciesInfo_atomFillings.cpp
	void estimateAtomEigs(); //!< If not read from file, estimate atomic eigenvalues from orbitals.
	void getAtom_nRadial(int spin, double magneticMoment, RadialFunctionG& nRadial) const; //!< Compute the atomic density per spin channel, given the magnetic moment
	
	friend class CommandIonSpecies;
	friend class CommandAddU;
	friend class CommandChargeball;
	friend class CommandWavefunction;
};


static EnumStringMap<SpeciesInfo::Constraint::ConstraintType> constraintTypeMap
(	SpeciesInfo::Constraint::None, "None",
	SpeciesInfo::Constraint::Linear, "Linear",
	SpeciesInfo::Constraint::Planar, "Planar"
);

#endif // JDFTX_ELECTRONIC_SPECIESINFO_H
