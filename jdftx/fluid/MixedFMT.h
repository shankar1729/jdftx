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

#ifndef JDFTX_FLUID_MIXEDFMT_H
#define JDFTX_FLUID_MIXEDFMT_H

#include <core/Operators.h>

//!@file MixedFMT.h
//!@brief Sphere mixture functional via (optionally soft) Fundamental Measure Theory


//! Returns the `White-Bear mark II' mixed sphere free energy/T given the weighted densities n*
//! and accumulates the gradients in grad_n*. Note that n1v and n2m are scalar
//! weighted densities, from which the vector and tensor weighted densities are obtained
//! internally by a gradient and traceless tensor second derivative respectively.
//! n2v is obtained as the negative gradient of n3. This is why n3, n1v and n2m
//! are passed in reciprocal space: they need fourier space processing for gradients etc.
double PhiFMT(const ScalarField& n0, const ScalarField& n1, const ScalarField& n2,
	const ScalarFieldTilde& n3tilde, const ScalarFieldTilde& n1vTilde, const ScalarFieldTilde& n2mTilde,
	ScalarField& grad_n0, ScalarField& grad_n1, ScalarField& grad_n2,
	ScalarFieldTilde& grad_n3tilde, ScalarFieldTilde& grad_n1vTilde, ScalarFieldTilde& grad_n2mTilde);

//! Returns the free energy density/T and accumulates derivatives
//! corresponding to PhiFMT() for the uniform fluid
double phiFMTuniform(double n0, double n1, double n2, double n3,
	double& grad_n0, double& grad_n1, double& grad_n2, double& grad_n3);


//! Bonding correction for tangentially bonded hard spheres
//! Rhm = Ra Rb /(Ra+Rb) is the harmonic sum of the sphere radii
//! scale is a scale factor for the correction (ratio of bond multiplicity to number of hard sphere sites in molecule)
//! n0mol is the suitably weighted partial measure-0 weighted density of this molecule
//! n2 and n3 are the usual FMT weighted densities.
//! Returns the free energy/T of bonding and accumulates gradients in grad_n*
//! Note that n3 is in fourier space for faster computation of n2v = -gradient n3
double PhiBond(double Rhm, double scale, const ScalarField& n0mol, const ScalarField& n2, const ScalarFieldTilde& n3tilde,
	ScalarField& grad_n0mol, ScalarField& grad_n2, ScalarFieldTilde& grad_n3tilde);

//! Returns the free energy density/T and accumulates derivatives
//! corresponding to PhiBond() for the uniform fluid
double phiBondUniform(double Rhm, double scale, double n0mol, double n2, double n3,
	double& grad_n0mol, double& grad_n2, double& grad_n3);

#endif // JDFTX_FLUID_MIXEDFMT_H
