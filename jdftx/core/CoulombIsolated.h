/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_COULOMBISOLATED_H
#define JDFTX_CORE_COULOMBISOLATED_H

#include <core/Coulomb.h>
#include <core/WignerSeitz.h>

//! @addtogroup LongRange
//! @{
//! @file CoulombIsolated.h Coulomb-interactions in 0D

//! Coulomb interaction for an isolated system (no periodicity), truncated on the Wigner-Seitz cell
class CoulombIsolated : public Coulomb
{
public:
	CoulombIsolated(const GridInfo& gInfoOrig, const CoulombParams& params);
protected:
	ScalarFieldTilde apply(ScalarFieldTilde&&) const;
	std::shared_ptr<Ewald> createEwald(matrix3<> R, size_t nAtoms) const;
	matrix3<> getLatticeGradient(const ScalarFieldTilde& X, const ScalarFieldTilde& Y) const;
private:
	WignerSeitz ws;
	RealKernel Vc;
	ManagedArray<symmetricMatrix3<>> Vc_RRT; //optional lattice derivative of kernel
	friend struct ExchangeEval;
};

//! Coulomb interaction for an isolated system (no periodicity), truncated on a sphere
class CoulombSpherical : public Coulomb
{
public:
	CoulombSpherical(const GridInfo& gInfoOrig, const CoulombParams& params);
protected:
	ScalarFieldTilde apply(ScalarFieldTilde&&) const;
	std::shared_ptr<Ewald> createEwald(matrix3<> R, size_t nAtoms) const;
	matrix3<> getLatticeGradient(const ScalarFieldTilde& X, const ScalarFieldTilde& Y) const;
private:
	WignerSeitz ws;
	double Rc; //!< radius of truncation (set to Wigner-Seitz in-radius if params.Rc=0)
	friend struct ExchangeEval;
	friend class BGW;
};

//! @}
#endif // JDFTX_CORE_COULOMBISOLATED_H
