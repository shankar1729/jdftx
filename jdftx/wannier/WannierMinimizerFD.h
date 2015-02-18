/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_WANNIER_WANNIERMINIMIZERFD_H
#define JDFTX_WANNIER_WANNIERMINIMIZERFD_H

#include <wannier/WannierMinimizer.h>

class WannierMinimizerFD : public WannierMinimizer
{
public:
	WannierMinimizerFD(const Everything& e, const Wannier& wannier);

	void initialize(int iSpin);
	double getOmega(bool grad);
	double getOmegaI(bool grad);
	WannierGradient precondition(const WannierGradient& grad);
	
	//!An edge of the k-mesh involved in the finite difference formula
	struct Edge
	{	double wb; //!< weight of neighbour
		vector3<> b; //!< displacement to neighbour
		unsigned ik; //!< index of neighbour in kMesh
		Kpoint point; //!< description of neighbour (source state, rotation, translation etc.)
		matrix M0; //!< initial overlap matrix for this pair
	};
	std::vector< std::vector<Edge> > edges;
	matrix kHelmholtzInv;
};

#endif //JDFTX_WANNIER_WANNIERMINIMIZERFD_H