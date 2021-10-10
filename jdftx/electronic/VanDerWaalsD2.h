/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Kendra Letchworth Weaver

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

#ifndef JDFTX_ELECTRONIC_VANDERWAALSD2_H
#define JDFTX_ELECTRONIC_VANDERWAALSD2_H

#include <electronic/VanDerWaals.h>
#include <core/RadialFunction.h>

//! @addtogroup LongRange
//! @{

//! DFT-D2 pair potential dispersion correction \cite Dispersion-D2
class VanDerWaalsD2 : public VanDerWaals
{
public:
	VanDerWaalsD2(const Everything &e, string reason="");
	~VanDerWaalsD2();
	
	const static int unitParticle = -1; //!< special atomic number used by some fluids: point particle with C6=1 J-nm^6/mol and R0=0
	
	//Implement virtual functions of VanDerWaals abstract base class
	virtual double getScaleFactor(string exCorrName, double scaleOverride=0.) const;
	virtual double energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT=0) const;

	//! Van der Waal correction to the interaction energy between the explicit atoms
	//! (from IonInfo) and the continuous fields Ntilde with specified atomic numbers.
	//! The gradient w.r.t site densities is accumulated to grad_Ntilde (if non-null),
	//! the negative gradient w.r.t discrete atom positions is accumulated to forces (if non-null)
	//! and the lattice gradient is accumulated to E_RRT (if non-null)
	double energyAndGrad(const std::vector< std::vector< vector3<> > >& atpos,
		const ScalarFieldTildeArray& Ntilde, const std::vector<int>& atomicNumber,
		const double scaleFac, ScalarFieldTildeArray* grad_Ntilde=0, struct IonicGradient* forces=0, matrix3<>* E_RRT=0) const;
	
	//! C6 and R0 parameters for the VDW interactions
	struct AtomParams
	{	double C6;
		double R0;
		//!Construct given C6 in J-nm^6/mol and R0 in Angstrom
		AtomParams(double SI_C6=0., double SI_R0=0.);
	};
	AtomParams getParams(int atomicNumber, int sp) const; //!< retrieve vdW parameters for an atom
	
private:
	std::vector<AtomParams> atomParams; //!< List of C6 coeficients and radii R0 for all atoms
	std::map<string,double> scalingFactor; //!< ExCorr dependent scale factor
	
	//! Get cached RadialFunctionG for interaction kernel between species of
	//! atomic numnbers Z1 and Z2. The radial function will be created (and cached)
	//! the first time this function is called with a specific pair of Z's.
	//! Note that both Z1 and Z2 should be supported non-zero atomic numbers (or the special value unitParticle)
	//! If sp is non-negative, then the corresponding species will be checked for overriding parameters
	const RadialFunctionG& getRadialFunction(int Z1, int Z2, int sp1, int sp2) const;
	
	std::map<std::pair<int,int>,RadialFunctionG> radialFunctions;
};

//! @}
#endif // JDFTX_ELECTRONIC_VANDERWAALSD2_H
