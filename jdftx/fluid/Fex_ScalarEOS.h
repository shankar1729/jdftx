/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#ifndef JDFTX_FLUID_FEX_SCALAREOS_H
#define JDFTX_FLUID_FEX_SCALAREOS_H

#include <fluid/Fex.h>

//! @addtogroup ClassicalDFT
//! @{
//! @file Fex_ScalarEOS.h ScalarEOS functional and helper classes

//! Abstract base class for the equation of state evaluator for ScalarEOS functionals
struct ScalarEOS
{	double sigmaEOS;
	virtual double vdwRadius() const=0;
	virtual void evaluate(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const=0;
	#ifdef GPU_ENABLED
	virtual void evaluate_gpu(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const=0;
	#endif
	ScalarEOS(double sigmaEOS) : sigmaEOS(sigmaEOS) {}
	virtual ~ScalarEOS() {}
};

//! Scalar EOS functional from \cite RigidCDFT and its extension from \cite PolarizableCDFT
class Fex_ScalarEOS : public Fex
{
public:
	Fex_ScalarEOS(const FluidMixture*, const FluidComponent*, const ScalarEOS& eos);
    virtual ~Fex_ScalarEOS();
	
	double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const;
	double computeUniform(const double* N, double* Phi_N) const;

private:
	const ScalarEOS& eos; double Vhs;
	RadialFunctionG fex_LJatt;
};

//! Jefferey-Austin equation of state for water
struct JeffereyAustinEOS : public ScalarEOS
{	JeffereyAustinEOS(double T, double sigmaEOS);
	double vdwRadius() const;
	void evaluate(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const;
	#ifdef GPU_ENABLED
	void evaluate_gpu(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const;
	#endif
private:
	std::shared_ptr<struct JeffereyAustinEOS_eval> eval;
};

//! Tao-Mason equation of state for moderately polar liquids
struct TaoMasonEOS : public ScalarEOS
{	TaoMasonEOS(double T, double Tc, double Pc, double omega, double sigmaEOS);
	double vdwRadius() const;
	void evaluate(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const;
	#ifdef GPU_ENABLED
	void evaluate_gpu(size_t nData, const double* N, double* Aex, double* Aex_N, double Vhs) const;
	#endif
private:
	std::shared_ptr<struct TaoMasonEOS_eval> eval;
};

//! @}
#endif // JDFTX_FLUID_FEX_SCALAREOS_H
