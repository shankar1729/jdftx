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

#ifndef JDFTX_ELECTRONIC_EXACTEXCHANGE_H
#define JDFTX_ELECTRONIC_EXACTEXCHANGE_H

#include <core/ScalarField.h>

class Everything;
class ColumnBundle;

//! @addtogroup ExchangeCorrelation
//! @{

//! Exact-exchange calculator
class ExactExchange
{
public:
	ExactExchange(const Everything& e);
	~ExactExchange();
	
	//! Compute scaled exact exchange energy with scale aXX and range omega
	//! (and optionally accumulate gradients) given fillings and wavefunctions.
	//! If prepareHamiltonian has been called with the same omega already,
	//! then the ACE representation is used to compute the energy and HC instead,
	//! except when lattice gradient or RPA-mode are requested which require full computation.
	//! RPA mode additionally requires electronic eigenvalues E.
	double operator()(double aXX, double omega,
		const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
		std::vector<ColumnBundle>* HC=0, matrix3<>* EXX_RRT=0,
		bool rpaMode=false, const std::vector<diagMatrix>* Hsub_eigs=0) const;
	
	//! Initialize the ACE (Adiabatic Compression of Exchange) representation in preparation for applyHamiltonian
	void prepareHamiltonian(double omega, const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C);
	
	//! Apply Hamiltonian using ACE representation initialized previously, and return the exchange energy contribution from current q.
	//! Note that fillings Fq are only used for computing the energy, and do not impact the Hamiltonian which only depends on F used in prepareHamiltonian().
	//! HCq must be allocated (non-null) in order to collect the Hamiltonian contribution, else only energy is returned.
	double applyHamiltonian(double aXX, double omega, int q, const diagMatrix& Fq,  const ColumnBundle& Cq, ColumnBundle& HCq) const;

	//! Add Hamiltonian in plane-wave basis (used for dense diagonalization for BGW):
	void addHamiltonian(double aXX, double omega, int q, matrix& H,
		const std::vector<int>& iRowsMine, const std::vector<int>& iColsMine) const;

private:
	const Everything& e;
	class ExactExchangeEval* eval; //!< opaque pointer to an internal computation class
};

//! @}
#endif // JDFTX_ELECTRONIC_EXACTEXCHANGE_H
