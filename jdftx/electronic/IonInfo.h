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

#ifndef JDFTX_ELECTRONIC_IONINFO_H
#define JDFTX_ELECTRONIC_IONINFO_H

#include <electronic/SpeciesInfo.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/IonicGaussianPotential.h>
#include <electronic/Metadynamics.h>
#include <core/matrix.h>
#include <core/ScalarField.h>
#include <core/Thread.h>

//! @addtogroup IonicSystem
//! @{
//! @file IonInfo.h Class IonInfo and related definitions

//! Coordinate system for ion positions
enum CoordsType {CoordsLattice, CoordsCartesian}; 

//! Coordinate system for force output:
enum ForcesOutputCoords { ForcesCoordsPositions, ForcesCoordsLattice, ForcesCoordsCartesian, ForcesCoordsContravariant };
static EnumStringMap<ForcesOutputCoords> forcesOutputCoordsMap(
	ForcesCoordsPositions, "Positions",
	ForcesCoordsLattice, "Lattice",
	ForcesCoordsCartesian, "Cartesian",
	ForcesCoordsContravariant, "Contravariant");

//! Check method used for determining whether pseudopotential cores overlap
enum coreOverlapCheck { additive, vector, none };

//! Style of vdW correction
enum VDWstyle { VDW_D2, VDW_D3 };


//! Container class for ionic system: collection of species, each with several atoms
class IonInfo
{
public:
	std::vector< std::shared_ptr<SpeciesInfo> > species; //!< list of ionic species
	std::vector<string> pspFilenamePatterns; //!< list of wildcards for pseudopotential sets
	
	CoordsType coordsType; //!< coordinate system for ionic positions etc.
	ForcesOutputCoords forcesOutputCoords; //!< coordinate system to print forces in
	coreOverlapCheck coreOverlapCondition; //! Check method used for determining whether pseudopotential cores overlap
	bool vdWenable; //!< whether vdW pair-potential corrections are enabled
	VDWstyle vdWstyle; //!< style of vdW correction for ions in electronic system (not used for fluid coupling)
	double vdWscale; //!< If non-zero, override the default scale parameter
	double ljOverride; //!< If non-zero, replace electronic DFT with LJ pair potential with rCut=ljOverride (for testing geometry optimization and dynamics algorithms only)
	std::vector<IonicGaussianPotential> ionicGaussianPotentials; //!< External Gaussian potentials and forces on atoms
	std::shared_ptr<MetadynamicsBond> metadynamicsBond; //!< Optional metadynamics on a bond degree of freedom
	
	IonicGradient forces; //!< forces at current atomic positions in latice coordinates
	matrix3<> stress; //!< stresses at current lattice geometry in Eh/a0^3 (only calculated if optimizing lattice or dumping stress)
	bool computeStress; //!< flag to control whether ionicEnergyAndGrad() computes the stress tensor in addition to forces
	diagMatrix thermostat; //!< optional thermostat internal degrees of freedom used for IonicDynamics
	diagMatrix barostat; //!< optional barostat internal degrees of freedom used for IonicDynamics
	
	ScalarFieldTilde Vlocps; //!< Net local pseudopotential
	ScalarFieldTilde rhoIon; //!< Total ionic charge density (with width ionWidth, used for interactions with fluid)
	ScalarFieldTilde nChargeball; //!< Extra electron density around ionic cores to keep fluid out (DEPRECATED)
	ScalarField nCore; //!< Core electron density for partial (nonlinear) core correction
	ScalarField tauCore; //!< Model for the KE density of the core (TF+vW applied to nCore) (used by meta-GGAs)
	
	IonInfo();
	
	void setup(const Everything&);
	void printPositions(FILE*) const; 
	bool checkPositions() const; //!< check for overlapping atoms, return true if okay
	double getZtot() const; //!< get total Z of all species and atoms
	
	//! Update Vlocps, rhoIon, nChargeball, nCore and the energies dependent only on ionic positions
	void update(class Energies&); 

	//! Return the total (free) energy and update the ionic gradient in IonInfo::forces
	//! If IonInfo::computeStress = true, also update the lattice gradient in IonInfo::stress
	double ionicEnergyAndGrad();

	//! Return the non-local pseudopotential energy due to a single state.
	//! Optionally accumulate the corresponding electronic gradient in HCq and ionic gradient in forces
	double EnlAndGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const std::vector<matrix>& VdagCq, std::vector<matrix>& HVdagCq) const;
	
	//! Accumulate pseudopotential dependent contribution to the overlap in OCq
	void augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq, std::vector<matrix>* VdagCq=0) const;
	
	//Multi-stage density augmentation and gradient propagation (see corresponding functions in SpeciesInfo)
	void augmentDensityInit() const; //!< initialize density augmentation
	void augmentDensitySpherical(const QuantumNumber& qnum, const diagMatrix& Fq, const std::vector<matrix>& VdagCq, const std::vector<matrix>* dVdagCqL = 0, const std::vector<matrix>* dVdagCqR = 0) const; //!< calculate density augmentation in spherical functions. Parameters dVdagCqL and/or dVdagCqR are needed for perturbation theory.
	void augmentDensityGrid(ScalarFieldArray& n) const; //!< propagate from spherical functions to grid
	void augmentDensityGridGrad(const ScalarFieldArray& E_n, IonicGradient* forces=0, matrix3<>* Eaug_RRT=0) const; //!< propagate grid gradients to spherical functions
	void augmentDensitySphericalGrad(const QuantumNumber& qnum, const std::vector<matrix>& VdagCq, std::vector<matrix>& HVdagCq) const; //!< propagate spherical function gradients to wavefunctions
	
	void setE_nAug(const std::vector<matrix> E_nAug);
	const std::vector<matrix> getE_nAug();
	
	void project(const ColumnBundle& Cq, std::vector<matrix>& VdagCq, matrix* rotExisting=0) const; //Update pseudopotential projections (optionally retain non-zero ones with specified rotation)
	void projectGrad(const std::vector<matrix>& HVdagCq, const ColumnBundle& Cq, ColumnBundle& HCq) const; //Propagate projected gradient (HVdagCq) to full gradient (HCq)
	
	//! Compute U corrections (DFT+U in the simplified rotationally-invariant scheme [Dudarev et al, Phys. Rev. B 57, 1505])
	//rhoAtom is a flat array of atomic density matrices per U type, with index order (outer to inner): species, Uparam(n,l), spin, atom
	size_t rhoAtom_nMatrices() const; //!< number of matrices in rhoAtom array
	void rhoAtom_initZero(std::vector<matrix>& rhoAtom) const; //!< initialize matrices of appropriate size to zero
	void rhoAtom_calc(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<matrix>& rhoAtom) const; //!< compute atomic density matrices
	double rhoAtom_computeU(const std::vector<matrix>& rhoAtom, std::vector<matrix>& U_rhoAtom) const; //!< compute U energy and gradient w.r.t atomic density matrices
	void rhoAtom_grad(const ColumnBundle& Cq, const std::vector<matrix>& U_rhoAtom, ColumnBundle& HCq) const; //!< propagate U_rhoAtom to wavefunction gradient (per k-point to enable band structure)
	void rhoAtom_forces(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, const std::vector<matrix>& U_rhoAtom,
		IonicGradient& forces, matrix3<>* EU_RRT) const; //!< propagate U_rhoAtom to forces (and optionally stress)
	
	//! Expectation value of [r_iDir,H] = D_iDir + nonlocal corrections
	matrix rHcommutator(const ColumnBundle &Y, int iDir, const matrix& YdagHY) const; 

	//! Expectation value of [r_iDir,H] = D_iDir + nonlocal corrections between different wavefunctions.
	//! Here, H12 = Y1 ^ H(Y2) is needed for ultrasoft overlap-augmentation corrections.
	//! Return a result for each Cartesian direction.
	vector3<matrix> rHcommutator(const ColumnBundle &Y1, const ColumnBundle &Y2, const matrix& H12) const; 

	int nAtomicOrbitals() const; //!< Get total number of atomic orbitals
	ColumnBundle getAtomicOrbitals(int q, bool applyO, int extraCols=0) const; //!< Get all atomic orbitals of a given state number q, optionally with operator O pre-applied (with room for extra columns if specified)
	
	//! Compute all pair-potential terms in the energy, forces or lattice derivative (E_RRT) (electrostatic, and optionally vdW)
	void pairPotentialsAndGrad(class Energies* ener=0, IonicGradient* forces=0, matrix3<>* E_RRT=0) const;
	
	//! Method for determining ion charge width
	enum IonWidthMethod
	{	IonWidthEcut, //!< determine ion width from Ecut
		IonWidthFFTbox, //!< determine ion width from grid spacing
		IonWidthManual //!< manually specify the ion width
	}
	ionWidthMethod; //!< method for determining ion charge width
	double ionWidth; //!< width for gaussian representation of nuclei
	bool shouldPrintForceComponents;

private:
	const Everything* e;
	ScalarFieldTilde rhoIonBare; //rhoIon without ionWidth required for stress calculation
	
	//! Compute pulay contributions to energy and optionally stress
	double calcEpulay(matrix3<>* E_RRT=0) const;
};

//! @}
#endif // JDFTX_ELECTRONIC_IONINFO_H
