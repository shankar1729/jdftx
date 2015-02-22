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

#ifndef JDFTX_CORE_COULOMBWIRE_H
#define JDFTX_CORE_COULOMBWIRE_H

#include <core/Coulomb.h>
#include <core/WignerSeitz.h>

//! Coulomb interaction for a 1D periodic system, truncated on the transverse Wigner-Seitz cell
class CoulombWire : public Coulomb
{
public:
	CoulombWire(const GridInfo& gInfoOrig, const CoulombParams& params);
protected:
	ScalarFieldTilde apply(ScalarFieldTilde&&) const;
	std::shared_ptr<Ewald> createEwald(matrix3<> R, size_t nAtoms) const;
private:
	WignerSeitz ws;
	RealKernel Vc;
};

//! Coulomb interaction for a 1D periodic system, truncated on a cylinder
class CoulombCylindrical : public Coulomb
{
public:
	CoulombCylindrical(const GridInfo& gInfoOrig, const CoulombParams& params);
protected:
	ScalarFieldTilde apply(ScalarFieldTilde&&) const;
	std::shared_ptr<Ewald> createEwald(matrix3<> R, size_t nAtoms) const;
private:
	WignerSeitz ws;
	double Rc; //!< radius of truncation (set to Wigner-Seitz in-radius if params.Rc=0)
	RealKernel Vc;
	friend class ExchangeEval;
};

#endif // JDFTX_CORE_COULOMBWIRE_H
