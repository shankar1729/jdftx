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
#include <core/DataMultiplet.h>
#include <core/matrix3.h>

//! Convert a complex wavefunction to a real one with optimum phase choice
//! Store the phase statistics before conversion in meanPhase and sigmaPhase
//! and the relative rms imaginary part truncated during conversion in rmsImagErr
//! (Useful for getting real wavefunctions in gamma point only calculations or Wannier functions)
void removePhase(size_t N, complex* data, double& meanPhase, double& sigmaPhase, double& rmsImagErr);

//----------------- Scalar field operators ---------------------------

DataGptr D(const DataGptr&, int iDir); //!< compute the gradient in the iDir'th cartesian direction
DataGptr DD(const DataGptr&, int iDir, int jDir); //!< second derivative along iDir'th and jDir'th cartesian directions

//! Multiply complex scalar field by Block phase for wave-vector k (in reciprocal lattice coordinates)
void multiplyBlochPhase(complexDataRptr&, const vector3<>& k);

//! Resample scalar field by gathering from coordinates rotated by point group element
//! specified in mesh coordinates. Hermitian conjugate of pointGroupScatter
DataRptr pointGroupGather(const DataRptr&, const matrix3<int>& mMesh);
complexDataRptr pointGroupGather(const complexDataRptr&, const matrix3<int>& mMesh);

//! Resample scalar field by scattering to coordinates rotated by point group element
//! specified in mesh coordinates. Hermitian conjugate of pointGroupGather
DataRptr pointGroupScatter(const DataRptr&, const matrix3<int>& mMesh);
complexDataRptr pointGroupScatter(const complexDataRptr&, const matrix3<int>& mMesh);


//! Create a spherically symmetric real scalar field centered at lattice coordinates r0, given its radial fourier transform f
DataRptr radialFunction(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0); //calls DataGptr radialFunctionG

//! Create a spherically symmetric scalar G-space kernel given its radial form f 
void radialFunctionG(const RadialFunctionG& f, RealKernel& Kernel);  //calls DataGptr radialFunctionG

//! Create a spherically symmetric G-space scalar field centered at lattice coordinates r0 given its radial form f 
DataGptr radialFunctionG(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0);

//! Convolve a scalar field by a radial function (preserve input)
DataGptr operator*(const RadialFunctionG&, const DataGptr&);

//! Convolve a scalar field by a radial function (destructible input)
DataGptr operator*(const RadialFunctionG&, DataGptr&&);

//! Convolve a vector field by a radial function (preserve input)
DataGptrVec operator*(const RadialFunctionG&, const DataGptrVec&);

//! Convolve a vector field by a radial function (destructible input)
DataGptrVec operator*(const RadialFunctionG&, DataGptrVec&&);


//------------------------------ ColumnBundle operators ---------------------------------

//! Return Idag V .* I C (evaluated columnwise)
ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const DataRptr& V); 

//! Return projection (I-P)Y with P=O*C*C^
ColumnBundle Pbar(const ColumnBundle &C,const ColumnBundle &Y); 

ColumnBundle L(const ColumnBundle &Y); //!< Apply Laplacian
ColumnBundle Linv(const ColumnBundle &Y); //!< Apply Laplacian inverse
ColumnBundle O(const ColumnBundle &Y); //!< Apply overlap 
ColumnBundle D(const ColumnBundle &Y, int iDir); //!< Compute the cartesian gradient of a column bundle in direction# iDir

//! Apply inverse kinetic preconditioner inv((k+G)^2/2)
ColumnBundle precond_inv_kinetic(const ColumnBundle &Y, double KErollover); 

ColumnBundle translate(ColumnBundle&&, vector3<> dr); //!< translate a column-bundle by dr in lattice coordinates (destructible input)
ColumnBundle translate(const ColumnBundle&, vector3<> dr); //!< translate a column-bundle by dr in lattice coordinates (preserve input)

ColumnBundle switchBasis(const ColumnBundle&, const Basis&); //!< return wavefunction projected to a different basis

//------------------------------ ColumnBundle reductions ---------------------------------

//! Return trace(F*X^Y)
complex traceinner(const diagMatrix &F, const ColumnBundle &X,const ColumnBundle &Y);

//! Returns diag((I*X)*F*(I*X)^) (Compute density from an orthonormal wavefunction ColumnBundle with some fillings F)
DataRptr diagouterI(const diagMatrix &F,const ColumnBundle &X);

#endif // JDFTX_ELECTRONIC_OPERATORS_H
