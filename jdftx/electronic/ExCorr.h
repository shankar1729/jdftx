/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_EXCORR_H
#define JDFTX_ELECTRONIC_EXCORR_H

#include <electronic/common.h>
#include <core/DataCollection.h>

//! Types of exchange correlation functionals
enum ExCorrType
{
	ExCorrLDA_PZ, //!< Perdew-Zunger LDA
	ExCorrLDA_PW, //!< Perdew-Wang LDA (original version used in PW91)
	ExCorrLDA_PW_prec, //!< Perdew-Wang LDA with higher precision constants used in PBE
	ExCorrLDA_VWN, //!< Vosko-Wilk-Nusair LDA
	ExCorrLDA_Teter, //!< Teter93 LSDA functional
	ExCorrGGA_PBE, //!< PBE GGA functional
	ExCorrGGA_PBEsol, //!< PBE GGA functional reparametrized for solids
	ExCorrGGA_PW91, //!< PW91 GGA functional
	ExCorrMGGA_TPSS, //!< TPSS meta-GGA functional
	ExCorrMGGA_revTPSS, //!< revised meta-GGA functional
#ifdef LIBXC_ENABLED
	ExCorrLibXC,
#endif
	ExCorrHYB_PBE0, //!< PBE0 Hybrid GGA functional
	ExCorrHYB_HSE06, //!< HSE06 Screened Hybrid GGA functional
	ExCorrHF //!< Hartree-Fock
};

//! Types of kinetic energy functionals
enum KineticType
{
	KineticNone, //!< No kinetic energy (default)
	KineticTF, //!< Thomas-Fermi LDA kinetic energy
	KineticVW, //!< von Weisacker GGA kinetic energy
	KineticPW91 //!< PW91 GGA kinetic energy
};


class ExCorr
{
public:
	ExCorr(); 
	void setup(const Everything&); //!< Initialize
	string getName() const; //!< Get a description of the DFT functional
	
	//! Compute the exchange-correlation energy (and optionally gradient) for a (spin) density n
	//! Include kinetic energy if includeKinetic is set to true
	//! Orbital KE density tau must be provided if needsKEdensity() is true (for meta GGAs)
	//! and the corresponding gradient will be returned in Vtau if non-null
	//! For metaGGAs, Vtau should be non-null if Vxc is non-null
	double operator()(const DataRptrCollection& n, DataRptrCollection* Vxc=0, bool includeKinetic=false,
		const DataRptrCollection* tau=0, DataRptrCollection* Vtau=0) const;
	
	//! Compute the exchange-correlation energy (and optionally gradient) for a unpolarized density n
	//! Include kinetic energy if includeKinetic is set to true
	//! Orbital KE density tau must be provided if needsKEdensity() is true (for meta GGAs)
	//! and the corresponding gradient will be returned in Vtau if non-null.
	//! For metaGGAs, Vtau should be non-null if Vxc is non-null
	double operator()(const DataRptr& n, DataRptr* Vxc=0, bool includeKinetic=false,
		const DataRptr* tau=0, DataRptr* Vtau=0) const;

	double exxFactor() const; //!< retrieve the exact exchange scale factor (0 if no exact exchange)
	double exxRange() const; //!< range parameter (omega) for screened exchange (0 for long-range exchange)
	bool needsKEdensity() const; //!< whether orbital KE density is required as an input (for meta GGAs)
	
private:
	const Everything* e;
	ExCorrType exCorrType;
	KineticType kineticType;
	string xcName; // short name of the functional (set by command elec-ex-corr)
	friend class CommandElecExCorr;
	friend class CommandFluidExCorr;
	
	double exxScale; //scale factor for exact exchange
	double exxOmega; //Range parameter for exact exchange
	std::shared_ptr<class FunctionalList> functionals; //List of component functionals
#ifdef LIBXC_ENABLED
	int xcExchange, xcCorr, xcExcorr; //LibXC codes for various functional components (0 if unused)
#endif
};

#endif // JDFTX_ELECTRONIC_EXCORR_H
