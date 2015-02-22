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

#include <electronic/common.h>
#include <core/ScalarField.h>
#include <vector>

class ExactExchange
{
public:
	ExactExchange(const Everything& e);
	~ExactExchange();
	
	//! Compute scaled exact exchange energy with scale aXX and range omega
	//! (and optionally accumulate gradients) given fillings and wavefunctions
	double operator()(double aXX, double omega,
		const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
		std::vector<ColumnBundle>* HC = 0) const;
private:
	const Everything& e;
	class ExactExchangeEval* eval; //!< opaque pointer to an internal computation class
};

#endif // JDFTX_ELECTRONIC_EXACTEXCHANGE_H
