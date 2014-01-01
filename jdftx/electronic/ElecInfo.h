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

#include <electronic/common.h>
#include <core/vector3.h>

//! Spin polarization options
enum SpinType {SpinNone, SpinZ}; 


class QuantumNumber
{
public:
	vector3<> k; //!< k-point wave vector
	int spin;  //!< possible spin orientation. up=1, down=-1, none=0
	double weight; //!< state weight (= 1x or 2x k-point weight depending on spintype)

	QuantumNumber() : spin(0) {}
	int index() const { return spin<0 ? 1 : 0; } //!< return the appropriate index into electron (spin) density/potential arrays
};

class ElecInfo
{
public:
	int nBands, nStates; //!< Number of bands and total number of states
	int qStart, qStop; //!< Range of states handled by current process (= 0 and nStates for non-MPI jobs)
	bool isMine(int q) const { return q>=qStart && q<qStop; } //!< check if state index is local
	int whose(int q) const; //!< find out which process this state index belongs to
	int qStartOther(int iProc) const { return iProc ? qStopArr[iProc-1] : 0; } //!< find out qStart for another process
	int qStopOther(int iProc) const { return qStopArr[iProc]; } //!< find out qStop for another process
	
	SpinType spinType; //!< tells us what sort of spins we are using if any
	bool spinRestricted; //!< whether the calculation is spin restricted
	double nElectrons; //!< the number of electrons = Sum w Tr[F]
	std::vector<QuantumNumber> qnums; //!< k-points, spins and weights for each state
	
	enum FillingsUpdate
	{	ConstantFillings, //!< constant fillings (T=0)
		FermiFillingsMix, //!< mix fermi functions every mixInterval iterations
		FermiFillingsAux //!< fillings are a fermi function of the (auxilliary) subspace hamiltonian (recommended)
	}
	fillingsUpdate;
	
	double kT; //!< Temperature for Fermi distribution of fillings
	double mu; //!< If NaN, fix nElectrons, otherwise fix/target chemical potential to this
	
	int mixInterval; //!< we recalc. fillings every so many iterations
	bool subspaceRotation; //!< whether subspace variables are required (either rotation or aux hamiltonian)
	
	bool hasU; //! Flag to check whether the calculation has a DFT+U self-interaction correction

	std::vector<std::tuple<int, int, double>> customFillings;
	
	ElecInfo();
	void setup(const Everything &e, std::vector<diagMatrix>& F, Energies& ener); //!< setup bands and initial fillings
	void printFillings(FILE* fp) const;
	void mixFillings(std::vector<diagMatrix>& F, Energies& ener); //!< Fermi fillings with mixing / mu control
	void updateFillingsEnergies(const std::vector<diagMatrix>& F, Energies&) const; //!< Calculate fermi fillings Legendre multipliers (TS/muN)

	//Fermi function utilities:
	double fermi(double mu, double eps) const { return 0.5*(1.-tanh(betaBy2*(eps-mu))); } //!< fermi function
	double fermiPrime(double mu, double eps) const { return -0.5*betaBy2/pow(cosh(betaBy2*(eps-mu)), 2); } //!< derivative of fermi function
	diagMatrix fermi(double mu, const diagMatrix& eps) const; //!< elementwise fermi function
	diagMatrix fermiPrime(double mu, const diagMatrix& eps) const; //!< elementwise fermi function derivative
	
	//! Propagate matrix gradient w.r.t F to gradient w.r.t. eps (in the basis where fillings are diagonal)
	matrix fermiGrad(double mu, const diagMatrix& eps, const matrix& gradF) const;

	//! Compute number of electrons for a fermi distribution with specified eigenvalues
	double nElectronsFermi(double mu, const std::vector<diagMatrix>& eps) const; 
	
	//! Find the chemical potential for which the fermi distribution with specified eigenvalues adds up to nElectrons
	double findMu(const std::vector<diagMatrix>& eps, double nElectrons) const; 

	//!Find the best fit chemical potential (and optionally the density of states) given fillings and eigenvalues
	double fitMu(const std::vector<diagMatrix>& F, const std::vector<diagMatrix>& eps, double* dndmu=0) const;
	
	void kpointsPrint(bool printSpin=false) const; //!< Output k-points, weights and optionally spins
	void kpointPrint(int q, bool printSpin=false) const; //!< Output k-points, weights and optionally spins
	
	int findHOMO(int q) const; //! Returns the band index of the Highest Occupied Kohn-Sham Orbital

	//Parallel I/O utilities for diagMatrix/matrix array (one-per-kpoint, with nBands rows and columns unless overridden):
	void read(std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const;
	void read(std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const;
	void write(const std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const;
	void write(const std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const;

private:
	const Everything* e;
	double betaBy2; //!< initialized to 0.5/kT
	std::vector<int> qStopArr; //!< array of qStop's for all processes (used by whose() to find locate process managing a specific q)
	
	//Initial fillings:
	string initialFillingsFilename; //!< filename for initial fillings (zero-length if none)
	int nBandsOld; //!<number of bands in file being read
	std::vector<double> qNet; //!< net excess electrons in each spin channel
	friend class CommandElecInitialFillings;
	friend class CommandElecInitialCharge;
	friend class CommandInitialState;
	friend class ElecVars;
	
	//Fillings mix / mu-controller properties:
	double fillingMixFraction; //!< amount of new fillings mixed with old fillings
	double dnPrev, muMeasuredPrev; //!<The change in number of electrons and measured mu on the previous fillings update
	double Cmeasured, Cweight; //!<The "measured" capacitance of the system, and the weight of the measurement
	double dnMixFraction; //!<Scale the ideal stepsize in n by this factor
	friend class CommandElecFermiFillings;
	friend class CommandTargetMu;
	
	//k-points:
	vector3<int> kfold; //!< kpoint fold vector
	friend class CommandKpointFolding;
	friend class Everything;
	void kpointsFold(); //!< Fold k-points by kfold
	void kpointsReduce(); //!< Reduce folded k-points under symmetries
};

#endif // JDFTX_ELECTRONIC_ELECINFO_H
