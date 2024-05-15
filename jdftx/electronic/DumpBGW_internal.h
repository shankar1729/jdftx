/*-------------------------------------------------------------------
Copyright 2019 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_DUMPBGW_INTERNAL_H
#define JDFTX_ELECTRONIC_DUMPBGW_INTERNAL_H

//! @addtogroup Output
//! @{
//! @file DumpBGW_internal.h Implementation internals for BGW output module

#include <core/vector3.h>

//! Parameters for BGW output
struct BGWparams
{	int nBandsDense; //!< if non-zero, use a dense ScaLAPACK solver to calculate more bands
	int blockSize; //!< block size for ScaLAPACK diagonalization
	int clusterSize; //!< maximum eigenvalue cluster size to allocate extra ScaLAPACK workspace for
	int nBandsV; //!< if non-zero, number of bands for Vxc and Vxx output
	bool saveVxc; //!< whether to write exchange-correlation matrix elements
	bool saveVxx; //!< whether to write exact-exchange matrix elements
	bool rpaExx; //!< whether to compute RPA-consistent exact-exchange energy
	bool offDiagV; //!< whether to write off-diagonal matrix elements of Vxc and/or Vxx (default: false)

	double EcutChiFluid; //!< KE cutoff for fluid polarizability output (enabled if non-zero)
	bool elecOnly; //!< whether to only output electronic polarizability of fluid (default: true)
	vector3<> q0; //!< zero wavevector replacement used for polarizability output
	double freqReMax_eV; //!< maximum real frequency in eV
	double freqReStep_eV; //!< real frequency grid spacing in eV
	double freqBroaden_eV; //!< broadening (imaginary part) of real frequency grid in eV
	int freqNimag; //!< number of imaginary frequencies
	double freqPlasma; //!< plasma frequency in Hartrees used in GW imaginary frequency grid, set to zero for RPA frequency grid
	
	double Ecut_rALDA; //!< KE cutoff (in Eh) for rALDA output (enabled if non-zero)
	double kFcut_rALDA; //!< kF cutoff (in 1/a0) for rALDA regularization (enabled if non-zero)
	bool kernelSym_rALDA; //!< whether to use kernel symmetrization for rALDA (wavevector symmetrization if false, the default)
	
	BGWparams() : nBandsDense(0), blockSize(32), clusterSize(10), nBandsV(0),
		saveVxc(true), saveVxx(false), rpaExx(false), offDiagV(false),
		EcutChiFluid(0.), elecOnly(true),
		freqReMax_eV(30.), freqReStep_eV(1.), freqBroaden_eV(0.1),
		freqNimag(25), freqPlasma(1.),
		Ecut_rALDA(0.), kFcut_rALDA(0.), kernelSym_rALDA(false)
	{}
};


#ifdef HDF5_ENABLED //BGW output requires HDF5

#include <electronic/Everything.h>
#include <core/H5io.h>

//! Helper class for DumpBGW
//Methods implemented in DumpBGW.cpp except where indicated otherwise
class BGW
{
	const Everything& e;
	const BGWparams& bgwp;
	const GridInfo& gInfo;
	const ElecInfo& eInfo;
	const ElecVars& eVars;
	const int nSpins; //!< number of spin channels
	const int nSpinor; //!< number of spinor component
	const int nReducedKpts; //!< number of reduced k-points
	std::vector<int> nBasis; //!< number of basis functions per reduced k-point
	std::vector<int> nBasisPrev; //!< cumulative number of basis functions before reduced k-point
	int nBasisMax; //!< maximum number of basis functions for any reduced k-point
	std::vector<vector3<>> k; //!< k-points in BGW convention [0,1)
	std::vector<vector3<int>> kOffset; //!< k offsets to switch from JDFTx to BGW convention
	std::vector<double> wk; //!< k-point weights
	std::vector<vector3<int>> iGrho; //!< charge-density G-vector mesh
	
	int nBands; //!< eInfo.nBands, or overridden by nBandsDense
	int nBandsV; //!< same as nBands, unless overridden by non-zero BGWParams::nBandsV
	std::vector<diagMatrix> E, F; //!< eigenvalues and fillings
	std::vector<matrix> VxcSub; //!< exchange-correlation matrix elements
	std::vector<matrix> VxxSub; //!< exact exchange matrix elements
	std::vector<diagMatrix> VxcDiag, VxxDiag; //!< diagonal parts of the above
	
	hid_t openHDF5(string fname) const; //!< Open HDF5 file for collective access
	void writeHeaderMF(hid_t fid) const; //!< Write common HDF5 header specifying the mean-field claculation for BGW outputs
	void writeHeaderEps(hid_t gidHeader, bool write_q0, string mode,
		std::vector<vector3<>>& q, std::vector<complex>& freq,
		std::vector<std::vector<vector3<int>>>& iGarr,
		std::vector<int>& nBasis, int& nBasisMax) const; //!< Initialize header for polarizabilities, along with q and G-space quantities
	
	void writeV(std::vector<matrix>& Vsub, std::vector<diagMatrix>& Vdiag, string fname) const; //!< Common matrix elements I/O for writeVxc and writeVxx

public:
	BGW(const Everything& e, const BGWparams& bgwp);
	void writeWfn() const; //!< Write wavefunction file
	void writeVxc() const; //!< Write exchange-correlation matrix elements
	void writeVxx() const; //!< Write exact exchange matrix elements
	
	//Implemented in DumpBGW_dense.cpp
	void denseWriteWfn(hid_t gidWfns); //!< Solve and output wavefunctions using ScaLAPACK
	
	//Implemented in DumpBGW_fluid.cpp
	void writeChiFluid(bool write_q0) const; //!< Write fluid polarizability (for q0 or for q != q0 depending on write_q0)
	
	//Implemented in DumpBGW.cpp
	void write_rALDA(bool write_q0) const; //!< Write rALDA kernel (for q0 or for q != q0 depending on write_q0)
};

#endif //HDF5_ENABLED

//! @}
#endif //JDFTX_ELECTRONIC_DUMPBGW_INTERNAL_H
