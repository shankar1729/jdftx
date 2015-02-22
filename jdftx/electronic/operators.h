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

//! @file operators.h Operators specific to the electronic code

#include <electronic/common.h>
#include <electronic/RadialFunction.h>
#include <core/VectorField.h>
#include <core/ScalarFieldArray.h>
#include <core/matrix3.h>

//! Convert a complex wavefunction to a real one with optimum phase choice
//! Store the phase statistics before conversion in meanPhase and sigmaPhase
//! and the relative rms imaginary part truncated during conversion in rmsImagErr
//! (Useful for getting real wavefunctions in gamma point only calculations or Wannier functions)
void removePhase(size_t N, complex* data, double& meanPhase, double& sigmaPhase, double& rmsImagErr);

//----------------- Scalar field operators ---------------------------

ScalarFieldTilde D(const ScalarFieldTilde&, int iDir); //!< compute the gradient in the iDir'th cartesian direction
ScalarFieldTilde DD(const ScalarFieldTilde&, int iDir, int jDir); //!< second derivative along iDir'th and jDir'th cartesian directions
ScalarFieldTildeArray lGradient(const ScalarFieldTilde&, int l); //!< spherical tensor gradient of order l (2l+1 outputs, multiplied by Ylm(Ghat) (iG)^l)
ScalarFieldTilde lDivergence(const ScalarFieldTildeArray&, int l); //!< spherical tensor divergence of order l (2l+1 inputs, multiplied by Ylm(Ghat) (iG)^l, and summed)


//! Multiply complex scalar field by Block phase for wave-vector k (in reciprocal lattice coordinates)
void multiplyBlochPhase(complexScalarField&, const vector3<>& k);

//! Resample scalar field by gathering from coordinates rotated by point group element
//! specified in mesh coordinates. Hermitian conjugate of pointGroupScatter
ScalarField pointGroupGather(const ScalarField&, const matrix3<int>& mMesh);
complexScalarField pointGroupGather(const complexScalarField&, const matrix3<int>& mMesh);

//! Resample scalar field by scattering to coordinates rotated by point group element
//! specified in mesh coordinates. Hermitian conjugate of pointGroupGather
ScalarField pointGroupScatter(const ScalarField&, const matrix3<int>& mMesh);
complexScalarField pointGroupScatter(const complexScalarField&, const matrix3<int>& mMesh);


//! Create a spherically symmetric real scalar field centered at lattice coordinates r0, given its radial fourier transform f
ScalarField radialFunction(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0); //calls ScalarFieldTilde radialFunctionG

//! Create a spherically symmetric scalar G-space kernel given its radial form f 
void radialFunctionG(const RadialFunctionG& f, RealKernel& Kernel);  //calls ScalarFieldTilde radialFunctionG

//! Create a spherically symmetric G-space scalar field centered at lattice coordinates r0 given its radial form f 
ScalarFieldTilde radialFunctionG(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0);

//! Convolve a scalar field by a radial function (preserve input)
ScalarFieldTilde operator*(const RadialFunctionG&, const ScalarFieldTilde&);

//! Convolve a scalar field by a radial function (destructible input)
ScalarFieldTilde operator*(const RadialFunctionG&, ScalarFieldTilde&&);

//! Convolve a vector field by a radial function (preserve input)
VectorFieldTilde operator*(const RadialFunctionG&, const VectorFieldTilde&);

//! Convolve a vector field by a radial function (destructible input)
VectorFieldTilde operator*(const RadialFunctionG&, VectorFieldTilde&&);


//------------------------------ ColumnBundle operators ---------------------------------

//! Return Idag V .* I C (evaluated columnwise)
//! The handling of the spin structure of V parallels that of diagouterI, with V.size() taking the role of nDensities
ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const ScalarFieldArray& V);

ColumnBundle L(const ColumnBundle &Y); //!< Apply Laplacian
ColumnBundle Linv(const ColumnBundle &Y); //!< Apply Laplacian inverse
ColumnBundle O(const ColumnBundle &Y, std::vector<matrix>* VdagY=0); //!< Apply overlap (and optionally retrieve pseudopotential projections for later reuse)
ColumnBundle D(const ColumnBundle &Y, int iDir); //!< Compute the cartesian gradient of a column bundle in direction# iDir
ColumnBundle DD(const ColumnBundle &Y, int iDir, int jDir); //!< Compute second spatial derivative of a column bundle along directions# iDir, jDir

//! Apply inverse kinetic preconditioner inv((k+G)^2/2)
ColumnBundle precond_inv_kinetic(const ColumnBundle &Y, double KErollover); 

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

#endif // JDFTX_ELECTRONIC_OPERATORS_H
