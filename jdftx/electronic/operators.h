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

#ifndef JDFTX_ELECTRONIC_OPERATORS_H
#define JDFTX_ELECTRONIC_OPERATORS_H

#include <electronic/ColumnBundle.h>
#include <core/RadialFunction.h>
#include <core/VectorField.h>
#include <core/ScalarFieldArray.h>
#include <core/matrix3.h>

//! @addtogroup Operators
//! @{
//! @file operators.h Operators used extensively in the electronic sector

//------------------------------ ColumnBundle operators ---------------------------------

//! Return Idag V .* I C (evaluated columnwise)
//! The handling of the spin structure of V parallels that of diagouterI, with V.size() taking the role of nDensities
ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const ScalarFieldArray& V);

ColumnBundle L(const ColumnBundle &Y); //!< Apply Laplacian
ColumnBundle Linv(const ColumnBundle &Y); //!< Apply Laplacian inverse
ColumnBundle O(const ColumnBundle &Y, std::vector<matrix>* VdagY=0); //!< Apply overlap (and optionally retrieve pseudopotential projections for later reuse)
ColumnBundle D(const ColumnBundle &Y, int iDir); //!< Compute the cartesian gradient of a column bundle in direction# iDir
ColumnBundle DD(const ColumnBundle &Y, int iDir, int jDir); //!< Compute second spatial derivative of a column bundle along directions# iDir, jDir

//! Apply inverse kinetic preconditioner (Roughly inv((k+G)^2/2)) in-place
void precond_inv_kinetic(ColumnBundle &Y, double KErollover); 

diagMatrix diagDot(const ColumnBundle& X, const ColumnBundle& Y); //!< compute diag(X^Y) efficiently (avoid the off-diagonals)
void precond_inv_kinetic_band(ColumnBundle& Y, const diagMatrix& KEref); //!< In-place inverse kinetic preconditioner with band-by-band KE reference (Used by BandDavidson)

ColumnBundle translate(ColumnBundle&&, vector3<> dr); //!< translate a column-bundle by dr in lattice coordinates (destructible input)
ColumnBundle translate(const ColumnBundle&, vector3<> dr); //!< translate a column-bundle by dr in lattice coordinates (preserve input)
void translateColumns(ColumnBundle&, const vector3<>* dr); //!< translate each column of a column bundle by a different dr (in-place)

ColumnBundle switchBasis(const ColumnBundle&, const Basis&); //!< return wavefunction projected to a different basis

//------------------------------ ColumnBundle reductions ---------------------------------

//! Return trace(F*X^Y)
complex traceinner(const diagMatrix &F, const ColumnBundle &X,const ColumnBundle &Y);

//! Returns diag((I*X)*F*(I*X)^) (Compute density from an orthonormal wavefunction ColumnBundle with some fillings F).
//! nDensities is the number of scalar field in the output and controls how spin/spinors are handled:
//!    1: return total density regardless of spin / spinor nature
//!    2: return spin density in X.qnum->index()'th component of the output (valid for non-spinor X only)
//!    4: return spin density-matrix (valid for spinor X only)
//! If gInfoOut is specified, function ensures that the output is changed to that grid (in case tighter wfns grid is in use)
ScalarFieldArray diagouterI(const diagMatrix &F,const ColumnBundle &X, int nDensities, const GridInfo* gInfoOut=0);

//! @}
#endif // JDFTX_ELECTRONIC_OPERATORS_H
