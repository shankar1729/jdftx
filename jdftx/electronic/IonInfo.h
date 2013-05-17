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

#include <electronic/common.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/matrix.h>
#include <core/Data.h>
#include <vector>
#include <core/Thread.h>

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

class IonInfo
{
public:
	std::vector< std::shared_ptr<SpeciesInfo> > species; //!< list of ionic species
	CoordsType coordsType; //!< coordinate system for ionic positions etc.
	ForcesOutputCoords forcesOutputCoords; //!< coordinate system to print forces in
	coreOverlapCheck coreOverlapCondition; //! Check method used for determining whether pseudopotential cores overlap
	bool vdWenable; //!< whether vdW pair-potential corrections are enabled
	double vdWscale; //!< If non-zero, override the default scale parameter
	
	IonicGradient forces; //!< forces at current atomic positions
	
	DataGptr Vlocps; //!< Net local pseudopotential
	DataGptr rhoIon; //!< Total ionic charge density (with width ionWidth, used for interactions with fluid)
	DataGptr nChargeball; //!< Extra electron density around ionic cores to keep fluid out (DEPRECATED)
	DataRptr nCore; //!< Core electron density for partial (nonlinear) core correction
	DataRptr tauCore; //!< Model for the KE density of the core (TF+vW applied to nCore) (used by meta-GGAs)
	
	IonInfo();
	
	void setup(const Everything&);
	void printPositions(FILE*) const; 
	bool checkPositions() const; //!< check for overlapping atoms, return true if okay
	double getZtot() const; //!< get total Z of all species and atoms
	double ionWidthMuCorrection() const; //!< correction to electron chemical potential due to finite ion width in fluid interaction
	
	//! Update Vlocps, rhoIon, nChargeball, nCore and the energies dependent only on ionic positions
	void update(Energies&); 

	//! Return the total (free) energy and calculate the ionic gradient (forces)
	double ionicEnergyAndGrad(IonicGradient& forces) const;

	//! Return the non-local pseudopotential energy due to a single state.
	//! Optionally accumulate the corresponding electronic gradient in HCq and ionic gradient in forces
	double EnlAndGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, IonicGradient* forces=0) const;
	
	//! Accumulate pseudopotential dependent contribution to the overlap in OCq
	void augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const;
	
	//! Accumulate the pseudopotential dependent contribution to the density in n
	void augmentDensity(const diagMatrix& Fq, const ColumnBundle& Cq, DataRptr& n) const;
	
	//! Propagate the gradient w.r.t n (Vscloc) to the gradient w.r.t Cq in HCq (if non-null)
	//! and to the ionic position gradient in forces (if non-null)
	//! The gradient w.r.t the overlap is also propagated to forces (gradCdagOCq must be non-null if forces is non-null)
	void augmentDensityGrad(const diagMatrix& Fq, const ColumnBundle& Cq, const DataRptr& Vscloc,
		ColumnBundle& HCq, IonicGradient* forces=0, const matrix& gradCdagOCq=matrix()) const;
	
	//! Compute U corrections (DFT+U in the simplified rotationally-invariant scheme [Dudarev et al, Phys. Rev. B 57, 1505])
	//! Also accumulate orbital gradients in HC, if non-null
	double computeU(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
		std::vector<ColumnBundle>* HC = 0, IonicGradient* forces=0) const;
	
	enum IonWidthMethod
	{	IonWidthEcut, //!< determine ion width from Ecut
		IonWidthFFTbox, //!< determine ion width from grid spacing
		IonWidthManual //!< manually specify the ion width
	}
	ionWidthMethod; //!< method for determining ion charge width
	double ionWidth; //!< width for gaussian representation of nuclei

private:
	const Everything* e;
	bool shouldPrintForceComponents;
	
	friend class CommandDebug;
	
	//! Compute all pair-potential terms in the energy or forces (electrostatic, and optionally vdW)
	void pairPotentialsAndGrad(Energies* ener=0, IonicGradient* forces=0) const;
};

#endif // JDFTX_ELECTRONIC_IONINFO_H
