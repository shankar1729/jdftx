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


class IonInfo
{
public:
	std::vector< std::shared_ptr<SpeciesInfo> > species; //!< list of ionic species
	CoordsType coordsType; //!< coordinate system for ionic positions etc.
	ForcesOutputCoords forcesOutputCoords; //!< coordinate system to print forces in
	
	IonicGradient forces; //!< forces at current atomic positions
	
	DataGptr Vlocps; //!< Net local pseudopotential
	DataGptr rhoIon; //!< Total ionic charge density (with width ionChargeWidth, used for interactions with fluid)
	DataGptr nChargeball; //!< Extra electron density around ionic cores to keep fluid out (DEPRECATED)
	DataRptr nCore; //!< Core electron density for partial (nonlinear) core correction
	DataRptr tauCore; //!< Model for the KE density of the core (TF+vW applied to nCore) (used by meta-GGAs)
	
	IonInfo();

	void setup(const Everything&);
	void printPositions(FILE*) const; 
	void checkPositions() const; //!< check for overlapping atoms
	double ElocCorrection() const; //!< get constant correction to Eloc due to nuclear width
	
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
	
	double GmaxNL; //!< maximum G extent for non-local projetcors (corresponds to Ecut)
	double GmaxLoc; //!< maximum G extent for local functions (corresponds to furthest fft-box vertex)

private:
	const Everything* e;
	bool shouldPrintForceComponents;
	
	enum IonWidthMethod
	{	IonWidthEcut, //!< determine ion width from Ecut
		IonWidthFFTbox, //!< determine ion width from grid spacing
		IonWidthManual //!< manually specify the ion width
	}
	ionWidthMethod; //!< method for determining ion charge width
	double ionChargeWidth; //!< width for gaussian representation of nuclei
	
	friend class IonicMinimizer;
	friend class SpeciesInfo;
	friend class CommandIonWidth;
	friend class CommandDebug;
	
	double sigma; //!< gaussian width for Ewald sums
	vector3<int> Nreal; //!< max unit cell indices for real-space part of Ewald sum
	vector3<int> Nrecip; //!< max unit cell indices for reciprocal-space part of Ewald sum
	
	//! Compute the ewald energy (nuclear-nuclear eletcrostatics) and optionally its gradient
	double ewaldAndGrad(IonicGradient* forces=0) const;
};

#endif // JDFTX_ELECTRONIC_IONINFO_H
