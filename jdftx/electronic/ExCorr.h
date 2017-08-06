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

#include <core/ScalarFieldArray.h>

//! @addtogroup ExchangeCorrelation
//! @{
//! @file ExCorr.h Class ExCorr and helpers

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
	ExCorrORB_GLLBsc, //!< GLLB-sc orbital-dependent potential (no total energy)
	ExCorrPOT_LB94, //!< Leeuwen-Baerends potentialfunctional (no total energy)
	ExCorrHYB_PBE0, //!< PBE0 Hybrid GGA functional
	ExCorrHYB_HSE06, //!< HSE06 Screened Hybrid GGA functional
	ExCorrHYB_HSE12, //! Reparametrized screened exchange functional for accuracy
	ExCorrHYB_HSE12s, //! Reparametrized screened exchange functional for minimum screening length
	ExCorrHF //!< Hartree-Fock
};

//! Types of kinetic energy functionals
enum KineticType
{
#ifdef LIBXC_ENABLED
	KineticLibXC,
#endif
	KineticNone, //!< No kinetic energy (default)
	KineticTF, //!< Thomas-Fermi LDA kinetic energy
	KineticVW, //!< von Weisacker GGA kinetic energy
	KineticPW91 //!< PW91 GGA kinetic energy
};

//! Which components to include in the results of ExCorr::operator()
struct IncludeTXC
{	bool T; //!< kinetic
	bool X; //!< exchange
	bool C; //!< correlation
	
	IncludeTXC(bool T=false, bool X=true, bool C=true) : T(T), X(X), C(C) {} //!< defaults to exchange-correlation without kinetic
};

//! Exchange-Correlation energy calculator
class ExCorr
{
public:
	ExCorr(ExCorrType exCorrType=ExCorrGGA_PBE, KineticType kineticType=KineticNone); 
	void setup(const Everything&); //!< Initialize
	string getName() const; //!< Get a description of the DFT functional
	
	//! Compute the exchange-correlation energy (and optionally gradient) for a (spin) density n.
	//! includeTXC selects which components to include in result (XC without kinetic by default).
	//! Orbital KE density tau must be provided if needsKEdensity() is true (for meta GGAs)
	//! and the corresponding gradient will be returned in Vtau if non-null
	//! For metaGGAs, Vtau should be non-null if Vxc is non-null
	double operator()(const ScalarFieldArray& n, ScalarFieldArray* Vxc=0, IncludeTXC includeTXC=IncludeTXC(),
		const ScalarFieldArray* tau=0, ScalarFieldArray* Vtau=0) const;
	
	//! Compute the exchange-correlation energy (and optionally gradient) for a unpolarized density n
	//! includeTXC selects which components to include in result (XC without kinetic by default).
	//! Orbital KE density tau must be provided if needsKEdensity() is true (for meta GGAs)
	//! and the corresponding gradient will be returned in Vtau if non-null.
	//! For metaGGAs, Vtau should be non-null if Vxc is non-null
	double operator()(const ScalarField& n, ScalarField* Vxc=0, IncludeTXC includeTXC=IncludeTXC(),
		const ScalarField* tau=0, ScalarField* Vtau=0) const;

	double exxFactor() const; //!< retrieve the exact exchange scale factor (0 if no exact exchange)
	double exxRange() const; //!< range parameter (omega) for screened exchange (0 for long-range exchange)
	bool needsKEdensity() const; //!< whether orbital KE density is required as an input (for meta GGAs)
	bool hasEnergy() const; //!< whether functional supports a total energy (if not, only usable in SCF, and no forces)
	
	//!Compute second derivatives of energy density w.r.t n and sigma=|grad n|^2 
	//! by finite difference (supported only for spin-unpolarized internal LDAs and GGAs).
	//! All sigma derivatives will be null on output for LDAs.
	//! The gradients will be set to zero for regions with n < nCut (useful to reduce numerical sensitivity in systems with empty space)
	void getSecondDerivatives(const ScalarField& n, ScalarField& e_nn, ScalarField& e_sigma, ScalarField& e_nsigma, ScalarField& e_sigmasigma, double nCut=1e-4) const;
	
	//! Abstract base class (interface specification) for orbital-dependent potential functionals
	struct OrbitalDep
	{	OrbitalDep(const Everything& e) : e(e) {}
		virtual ~OrbitalDep() {}
		virtual bool ignore_nCore() const=0; //!< Whether partial cores need to be ignored for this functional
		virtual ScalarFieldArray getPotential() const=0; //!< Return orbital-dependent portion of potential (obtains any necessary electronic property directly from ElecVars / ElecInfo)
		virtual void dump() const=0; //!< Dump any functional-specific quantities
	protected:
		const Everything& e;
	};
	std::shared_ptr<OrbitalDep> orbitalDep; //optional orbital-dependent potential functional
	
private:
	const Everything* e;
	ExCorrType exCorrType;
	KineticType kineticType;
	string xcName; // short name of the functional (set by command elec-ex-corr)
	friend struct CommandElecExCorr;
	friend struct CommandFluidExCorr;
	
	double exxScale; //scale factor for exact exchange
	double exxOmega; //Range parameter for exact exchange
	
	double exxScaleOverride, exxOmegaOverride; //override default values of EXX scale and omega
	friend struct CommandExchangeParameters; //command that sets the override parameters
	
	std::shared_ptr<struct FunctionalList> functionals; //List of component functionals
#ifdef LIBXC_ENABLED
	int xcExchange, xcCorr, xcExcorr, xcKinetic; //LibXC codes for various functional components (0 if unused)
#endif
};

//! @}
#endif // JDFTX_ELECTRONIC_EXCORR_H
