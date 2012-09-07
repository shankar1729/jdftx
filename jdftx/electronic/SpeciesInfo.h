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

#ifndef JDFTX_ELECTRONIC_SPECIESINFO_H
#define JDFTX_ELECTRONIC_SPECIESINFO_H

#include <electronic/common.h>
#include <electronic/RadialFunction.h>
#include <electronic/matrix.h>
#include <core/Data.h>
#include <core/vector3.h>
#include <core/string.h>




class SpeciesInfo
{
public:

	double Z; //!< Nuclear charge of the species
	string name; //!< Identifier

	std::vector<vector3<> > atpos; //!< array of atomic positions of this species
	#ifdef GPU_ENABLED
	vector3<> *atposGpu; //!< copy of atomic positions on the gpu
	void sync_atposGpu(); //!< update changes in atpos to atposGpu (call whenever atpos is changed)
	#endif
	vector3<>* atposPref; //!< points to atposGpu in GPU mode and atpos otherwise
	
	double dE_dnG; //!< Derivative of [total energy per atom] w.r.t [nPlanewaves per unit volume] (for Pulay corrections)
	double mass; //!< ionic mass (currently unused)	

	struct Constraint
	{	
		double moveScale;
		vector3<> x; //! The direction or plane normal of the constraint
		enum ConstraintType {None, Linear, Planar} type; //! Type of the constraint that is being applied to the ion
		
		vector3<> operator()(const vector3<>& grad);
		bool isEquivalent(const Constraint& otherConstraint, const matrix3<>& transform) const;
		int getDimension() const;
		void print(FILE* fp) const;
	};
	std::vector<Constraint> constraints; //! Constraint that is on the motion of the ions
	
	SpeciesInfo();
	~SpeciesInfo();
	void setup(const Everything&);
	void print(FILE* fp) const;

	//! Return non-local energy for this species and quantum number q an doptionally accumulate
	//! electronic gradient in HCq (if non-null) and ionic gradient in forces (if non-null)
	double EnlAndGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, std::vector<vector3<> >* forces) const;
	
	//! Accumulate pseudopotential contribution to the overlap in OCq
	void augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const;
	
	//! Accumulate the pseudopotential dependent contribution to the density in n
	void augmentDensity(const diagMatrix& Fq, const ColumnBundle& Cq, DataRptr& n) const;
	
	//! Propagate the gradient w.r.t n (Vscloc) to the gradient w.r.t Cq in HCq (if non-null)
	//! and to the ionic position gradient in forces (if non-null)
	//! The gradient w.r.t the overlap is also propagated to forces (gradCdagOCq must be non-null if forces is non-null)
	void augmentDensityGrad(const diagMatrix& Fq, const ColumnBundle& Cq, const DataRptr& Vscloc,
		ColumnBundle& HCq, std::vector<vector3<> >* forces=0, const matrix& gradCdagOCq=matrix()) const;
	
	//! Calculate atomic orbitals (store in Y with an optional column offset)
	void setAtomicOrbitals(ColumnBundle& Y, int colOffset=0) const;
	//! Store a single atomic orbital (iAtom'th atom, n'th shell of angular momentum l with specified m value) in col'th column of Y:
	void setAtomicOrbital(ColumnBundle& Y, int col, unsigned iAtom, unsigned n, int l, int m) const;
	int nAtomicOrbitals() const; //!< return number of atomic orbitals in this species (all atoms)
	int lMaxAtomicOrbitals() const; //!< return maximum angular momentum in available atomic orbitals
	int nAtomicOrbitals(int l) const; //!< return number of atomic orbitals of given l (per atom)
	int atomicOrbitalOffset(unsigned iAtom, unsigned n, int l, int m) const; //!< offset of specified atomic orbital in output of current species
	
	//! Binary-write the PAW projector matrices (if any) for a particular state, looped over atoms, l and then m
	void writeProjectors(const ColumnBundle& Cq, FILE* fp) const;
	
	//! Add contributions from this species to Vlocps, rhoIon, nChargeball and nCore/tauCore (if any)
	void updateLocal(DataGptr& Vlocps, DataGptr& rhoIon, DataGptr& nChargeball,
		DataGptr& nCore, DataGptr& tauCore) const; 
	
	//! Return the local forces (due to Vlocps, rhoIon, nChargeball and nCore/tauCore)
	std::vector< vector3<> > getLocalForces(const DataGptr& ccgrad_Vlocps, const DataGptr& ccgrad_rhoIon,
		const DataGptr& ccgrad_nChargeball, const DataGptr& ccgrad_nCore, const DataGptr& ccgrad_tauCore) const;
	
private:
	RadialFunctionG VlocRadial; //!< local pseudopotential
	RadialFunctionG nCoreRadial; //!< core density for partial core correction
	RadialFunctionG tauCoreRadial; //!< core KE density for partial core correction with meta-GGAs

	std::vector< std::vector<RadialFunctionG> > VnlRadial; //!< non-local projectors (outer index l, inner index projetcor)
	std::vector<matrix> Mnl; //!< nonlocal pseudopotential projector matrix (indexed by l)
	
	std::vector<matrix> Qint; //!< overlap augmentation matrix (indexed by l, empty if no augmentation)
	
	struct QijIndex
	{	int l1, p1; //!< Angular momentum and projector index for channel i
		int l2, p2; //!< Angular momentum and projector index for channel j
		int l; //!< net angular momentum (l varies from |l1-l2| to (l1+l2) in increments of 2)
		bool operator<(const QijIndex&) const; //!< comparison that ensures (l1,p1)<-->(l2,p2) symmetry
	private:
		void sortIndices(); //!< swap (l1,p1)<-->(l2,p2) indices if needed to bring to the upper triangular part
	};
	std::map<QijIndex,RadialFunctionG> Qradial; //!< radial functions for density augmentation
	
	bool readProjectors; //!< whether to read PAW projectors from Pot format files
	std::vector< std::vector<RadialFunctionG> > projRadial; //!< PAW projectors (outer index l, inner index projetcor)

	std::vector<std::vector<RadialFunctionG> > psiRadial; //!< radial part of the atomic orbitals (outer index l, inner index shell)
	
	enum PseudopotentialFormat
	{	Pot, //!< Old DFT++ format (.pot files)
		Cpi, //!< FHI98 format (.cpi files)
		Fhi, //!< FHI format with ABINIT header (.fhi files)
		Uspp //!< USPP format ultrasoft pseudopotential
	}
	pspFormat;
	string potfilename, pulayfilename;
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
	
	friend class CommandIonSpecies;
	friend class CommandChargeball;
	friend class CommandWavefunction;
	
};


static EnumStringMap<SpeciesInfo::Constraint::ConstraintType> constraintTypeMap
(	SpeciesInfo::Constraint::None, "None",
	SpeciesInfo::Constraint::Linear, "Linear",
	SpeciesInfo::Constraint::Planar, "Planar"
);

#endif // JDFTX_ELECTRONIC_SPECIESINFO_H
