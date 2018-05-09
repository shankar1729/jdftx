/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#ifndef JDFTX_ELECTRONIC_ELECINFO_H
#define JDFTX_ELECTRONIC_ELECINFO_H

#include <core/vector3.h>
#include <core/MPIUtil.h>

class matrix;
class diagMatrix;
class Energies;
class Everything;

//! @addtogroup ElecSystem
//! @{
//! @file ElecInfo.h

//! Spin polarization options
enum SpinType {
	SpinNone, //!< unpolarized
	SpinZ, //!< spin-polarized
	SpinVector, //!< noncollinear magnetism (supports spin-orbit)
	SpinOrbit //!< noncollinear but unpolarized (spin-orbit in nonmagnetic systems)
};

//! Quantum number of Bloch state: k-point and spin
class QuantumNumber
{
public:
	vector3<> k; //!< k-point wave vector
	int spin;  //!< possible spin orientation. up=1, down=-1, none=0
	double weight; //!< state weight (= 1x or 2x k-point weight depending on spintype)
	
	QuantumNumber() : spin(0) {}
	int index() const { return spin<0 ? 1 : 0; } //!< return the appropriate index into electron (spin) density/potential arrays
};

//! Conversion function needed for PeriodicLookup<QuantumNumber> used for k-point reduction
inline vector3<> getCoord(const QuantumNumber& qnum) { return qnum.k; }

//! Electronic system auxiliary information, fillings and related functions
class ElecInfo
{
public:
	int nBands, nStates; //!< Number of bands and total number of states
	int nDensities, spinWeight, qWeightSum; //!< number of density components, spin weight factor (= max occupation per state) and sum of k-point weights
	int qStart, qStop; //!< Range of states handled by current process (= 0 and nStates for non-MPI jobs)
	bool isMine(int q) const { return qDivision.isMine(q); } //!< check if state index is local
	int whose(int q) const { return qDivision.whose(q); } //!< find out which process this state index belongs to
	int qStartOther(int iProc) const { return qDivision.start(iProc); } //!< find out qStart for another process
	int qStopOther(int iProc) const { return qDivision.stop(iProc); } //!< find out qStop for another process
	
	SpinType spinType; //!< type of spin treatment
	double nElectrons; //!< the number of electrons = Sum w Tr[F]
	std::vector<QuantumNumber> qnums; //!< k-points, spins and weights for each state
	
	bool isNoncollinear() const { return spinType==SpinVector || spinType==SpinOrbit; }
	int spinorLength() const { return isNoncollinear() ? 2 : 1; }
	int nSpins() const { return spinType==SpinZ ? 2 : 1; }
	vector3<int> kFoldingCount() const { return kfold; }
	
	//! Fillings update mode
	enum FillingsUpdate
	{	FillingsConst, //!< constant fillings (T=0)
		FillingsHsub //!< fillings are a function of subspace Hamiltonian
	}
	fillingsUpdate; //!< fillings update mode
	bool scalarFillings; //!< whether fillings are scalar (equal for all bands) at all quantum numbers
	
	//! Smearing options
	enum SmearingType
	{	SmearingFermi, //!< Fermi-Dirac smearing
		SmearingGauss, //!< Gaussian smearing
		SmearingCold //!< Cold smearing
	}
	smearingType; //!< smearing option
	double smearingWidth; //!< Smearing width (temperature in the Fermi case)
	double mu; //!< If NaN, fix nElectrons, otherwise fix/target chemical potential to this
	double Bz; //!< If NaN, fix magnetization, otherwise fix/target magnetic field to this value
	bool muLoop; //!< Whether to optimize mu in an outer loop over fixed charge calculations
	
	bool hasU; //! Flag to check whether the calculation has a DFT+U self-interaction correction

	string initialFillingsFilename; //!< filename for initial fillings (zero-length if none)
	
	ElecInfo();
	void setup(const Everything &e, std::vector<diagMatrix>& F, Energies& ener); //!< setup bands and initial fillings
	void printFillings(FILE* fp) const;
	void smearReport(const double* muOverride=0) const; //Smearing report (compute mu from eigenvalues in eVars if muOverride not provided)
	void updateFillingsEnergies(const std::vector<diagMatrix>& eps, Energies&) const; //!< Calculate variable fillings Legendre multipliers (TS/muN)

	//Smearing function utilities:
	inline double muEff(double mu, double Bz, int q) const { return mu + Bz*qnums[q].spin; } //!< effective mu for each spin
	double smear(double mu, double eps) const; //!< smearing function
	double smearPrime(double mu, double eps) const; //!< derivative of smearing function
	double smearEntropy(double mu, double eps) const; //!< entropy associated with smearing function
	diagMatrix smear(double mu, const diagMatrix& eps) const; //!< elementwise smearing function
	diagMatrix smearPrime(double mu, const diagMatrix& eps) const; //!< elementwise smearing function derivative
	diagMatrix smearEntropy(double mu, const diagMatrix& eps) const; //!< element-wise smearing function entropy
	
	//! Propagate matrix gradient w.r.t F to gradient w.r.t. eps (in the basis where fillings are diagonal)
	matrix smearGrad(double mu, const diagMatrix& eps, const matrix& gradF) const;

	//! Compute number of electrons for the smearing function with specified eigenvalues
	//! If magnetization is constrained, bisect on corresponding Lagrange multiplier Bz, and retrieve it too
	double nElectronsCalc(double mu, const std::vector<diagMatrix>& eps, double& Bz) const; 
	
	//! Find the chemical potential for which the smearing function with specified eigenvalues adds up to nElectrons
	//! If magnetization is constrained, retrieve corresponding Lagrange multiplier Bz as well
	double findMu(const std::vector<diagMatrix>& eps, double nElectrons, double& Bz) const; 
	
	void kpointsPrint(FILE* fp, bool printSpin=false) const; //!< Output k-points, weights and optionally spins
	void kpointPrint(FILE* fp, int q, bool printSpin=false) const; //!< Output k-points, weights and optionally spins
	
	int findHOMO(int q) const; //!< Returns the band index of the Highest Occupied Kohn-Sham Orbital

	//Parallel I/O utilities for diagMatrix/matrix array (one-per-kpoint, with nBands rows and columns unless overridden):
	void read(std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const; //!< parallel read array of diagonal matrices
	void read(std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const; //!< parallel read array of matrices
	void write(const std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const; //!< parallel write array of diagonal matrices
	void write(const std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const;  //!< parallel write array of matrices

	//Parallel I/O utilities for ColumnBundle array (defined in COlumnBUndle.cpp):
	struct ColumnBundleReadConversion //!< Utility to convert columnbundle basis / bands
	{	bool realSpace; //!< whether to read realspace wavefunctions
		int nBandsOld; //!< nBands for the input wavefunction
		double Ecut, EcutOld; //!< Ecut for the current calculation and input wavefunction in fourier space
		vector3<int> S_old; //!< fftbox size for the input wavefunction in double space
		ColumnBundleReadConversion();
	};
	void read(std::vector<class ColumnBundle>&, const char *fname, const ColumnBundleReadConversion* conversion=0) const; //!< Read array of columnbundles, optionally with conversion
	void write(const std::vector<class ColumnBundle>&, const char *fname) const; //!< write an array of columnbundles to file

private:
	const Everything* e;
	TaskDivision qDivision; //!< MPI division of k-points
	
	//Initial fillings:
	int nBandsOld; //!<number of bands in file being read
	double Qinitial, Minitial; //!< net excess electrons and initial magnetization
	friend struct CommandElecInitialFillings;
	friend struct CommandElecInitialCharge;
	friend struct CommandElecInitialMagnetization;
	friend struct CommandInitialState;
	friend class ElecVars;
	friend struct LCAOminimizer;
	
	//!Calculate nElectrons and return magnetization at given mu, Bz and eigenvalues eps
	double magnetizationCalc(double mu, double Bz, const std::vector<diagMatrix>& eps, double& nElectrons) const; 
	
	//k-points:
	vector3<int> kfold; //!< kpoint fold vector
	friend struct CommandKpointFolding;
	friend class Everything;
	friend class Phonon;
	void kpointsFold(); //!< Fold k-points by kfold
	void kpointsReduce(); //!< Reduce folded k-points under symmetries
};

//! @}
#endif // JDFTX_ELECTRONIC_ELECINFO_H
