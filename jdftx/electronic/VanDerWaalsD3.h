/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman.

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

#ifndef JDFTX_ELECTRONIC_VANDERWAALSD3_H
#define JDFTX_ELECTRONIC_VANDERWAALSD3_H

#include <electronic/VanDerWaals.h>
#include <core/RadialFunction.h>

//! @addtogroup LongRange
//! @{


namespace D3
{	//! Parameter set per atom type
	struct AtomParams
	{	int Z; //!< atomic number
		double sqrtQ; //!< sqrt(Q) involved in computing C8 from C6
		double k2Rcov; //!< Covalent radius in bohrs, scaled by k2 = 4/3
		std::vector<double> CN; //!< reference coordination numbers
		std::vector<int> iCN; //!< index of each CN term in coefficient array
		
		matrix getL(double observedCN, matrix& Lprime) const; //!< compute weights and derivative of each reference CN at observed CN
	};
	
	//! Parameter set per pair of atom types
	struct PairParams
	{	double R0; //!< sum of cutoff radii for pair of atoms (in bohrs)
		matrix C6; //!< C6 coefficients at all pairs of reference CN's of the two atom types
	};
}


//! DFT-D3 pair potential dispersion correction \cite Dispersion-D3
class VanDerWaalsD3 : public VanDerWaals
{
public:
	VanDerWaalsD3(const Everything &e);
	
	//Implement virtual functions of VanDerWaals abstract base class
	virtual double getScaleFactor(string exCorrName, double scaleOverride=0.) const;
	virtual double energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRTptr=0) const;

private:
	double s6, s8; //scale factors for r^-6 and r^-8 terms corresponding to e.exCorr
	double sr6, sr8; //factors in damping of r^-6 and r^-8 terms corresponding to e.exCorr
	std::vector<D3::AtomParams> atomParams; //!< parameters per atom type
	std::vector<std::vector<D3::PairParams>> pairParams; //!< parameters per pair of atom types
	
	void computeCN(const std::vector<Atom>& atoms, std::vector<double>& CN) const; //!< compute coordination numbers
	void propagateCNgradient(const std::vector<Atom>& atoms, const std::vector<double>& E_CN,
		std::vector<vector3<>>& forces, matrix3<>* E_RRT) const; //!< propagate CN gradient to forces and stresses

	void report(const std::vector<double>& result, string name,
		const std::vector<Atom>& atoms, const char* fmt=" %.3f") const; //!<report per-atom quantity
	void setNeighborSampling(double rCut, vector3<int>& S, size_t& iStart, size_t& iStop) const; //!< get neighbor cells within rCut
};

//! @}
#endif // JDFTX_ELECTRONIC_VANDERWAALSD3_H
