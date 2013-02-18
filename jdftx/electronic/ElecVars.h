/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
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

#ifndef JDFTX_ELECTRONIC_ELECVARS_H
#define JDFTX_ELECTRONIC_ELECVARS_H

#include <electronic/common.h>
#include <electronic/FluidSolver.h>
#include <core/DataCollection.h>
#include <string>
#include <memory>

class ElecVars
{
public:
	//Independent variables:
	std::vector<ColumnBundle> Y; //!< unconstrained electronic wavefunctions
	std::vector<matrix> B; //!< subspace rotation / auxilliary hamiltonian
	double subspaceRotationFactor; //!< preconditioning factor for subspace rotations / aux hamiltonian relative to wavefunctions
	
	//Derived quantities:
	std::vector<matrix> B_evecs; //!<  eigenvectors of B[q] in columns
	std::vector<diagMatrix> B_eigs; //!< eigenvalues of B[q]

	std::vector<matrix> Hsub; //!< Subspace Hamiltonian:  Hsub[q]=C[q]^H*C[q]
	std::vector<matrix> Hsub_evecs; //!< eigenvectors of Hsub[q] in columns
	std::vector<diagMatrix> Hsub_eigs; //!< eigenvalues of Hsub[q]
	
	std::vector<ColumnBundle> C; //!< orthonormal wavefunctions (after appropriate subspace rotation)
	std::vector<diagMatrix> F;  //!< the fillings (diagonal matrices) for each state
	
	std::vector<matrix> grad_CdagOC; //!< gradient w.r.t overlap (required for forces when O is atom dependent)
	
	//Densities and potentials:
	DataRptrCollection n; //!< electron density (single DataR) or spin density (two DataR's)
	DataRptrCollection get_nXC() const; //!< return the total (spin) density including core contributions
	DataRptr get_nTot() const { return n.size()==1 ? n[0] : n[0]+n[1]; } //!< return the total electron density (even in spin polarized situations)
	
	DataRptrCollection tau; //!< kinetic energy density including tauCore, if present (computed if a meta GGA is being used)
	
	DataGptr d_vac; //!< electrostatic potential due to explicit electrons (and external charge, if any)
	DataGptr d_fluid; //!< electrostatic potential due to fluid
	DataGptr V_cavity; //!< non-electrostatic potential on electrons due to fluid
	IonicGradient fluidForces; //!< Extra forces due to fluid (besides contributions through d_fluid and V_cavity)
	
	DataRptrCollection Vscloc; //! Local part of (optionally spin-dependent) self-consistent potential
	DataRptrCollection Vtau; //! Gradient w.r.t kinetic energy density (if meta-GGA)
	
	//Fluid properties
	FluidType fluidType;
	FluidSolverParams fluidParams;
	std::shared_ptr<FluidSolver> fluidSolver;
	string fluidInitialStateFilename;

	//Wavefunction initialization:
	string wfnsFilename; //!< file to read wavefunctions from
	bool readWfnsRealspace; //!< whether to read realspace wavefunctions
	int nBandsOld; //!< nBands for the input wavefunction
	double EcutOld; BasisKdep kdepOld; //!< Ecut and basis-dependent flag  for the input wavefunction in fourier space
	int NxOld, NyOld, NzOld; //!< fftbox size for the input wavefunction in double space
	bool isRandom; //!< indicates whether the electronic state is random (not yet minimized)

	//Auxiliary hamiltonian initialization
	string HauxFilename; //!< file to read auxilliary hamiltonian (B) from (used only for FermiFillingsAux mode)
	bool HauxInitialized; //!< whether Haux has been read in/computed
	
	double overlapCondition; //!< Current condition number of the orbital overlap matrix (over all states)
	
	std::vector<string> nFilename; //!< file(s) to read electron (spin) density from

	ElecVars();
	void setup(const Everything &everything);

	//! Compute the terms written as a functional of the electronic density, and its gradient i.e. Vscloc
	void EdensityAndVscloc(Energies& ener);
	
	//! Update and return the electronic system/band structure energy.
	//! Optionally compute the gradient, preconditioned gradient and/or the subspace hamiltonian
	double elecEnergyAndGrad(Energies& ener, ElecGradient* grad=0, ElecGradient* Kgrad=0, bool calc_Hsub=false); 
	
	//! Set Y and C to eigenvectors of the subspace hamiltonian
	void setEigenvectors(); 
	
	//! Return the number of occupied bands (f > occupiedThrehsold) for a given state
	int nOccupiedBands(int q) const; 

	//!Compute the kinetic energy density
	DataRptrCollection KEdensity() const;
	
private:
	const Everything* e;
	
	std::vector<string> VexternalFilename; //!< external potential filename (read in real space)
	DataRptrCollection Vexternal; //!< external potential
	friend class CommandVexternal;
	friend class InverseKohnSham; //!< Adjusts Vexternal to produce target electron density
	
	string rhoExternalFilename; //!< external charge filename
	DataGptr rhoExternal; //!< external charge density
	friend class CommandRhoExternal;

	bool initLCAO; //!< initialize wave functions using linear combinations of atomic orbitals
	int LCAO(); //!< Initialize LCAO wavefunctions (returns the number of bands initialized)
	friend class CommandWavefunction;
	friend class Dump;
};
#endif // JDFTX_ELECTRONIC_ELECVARS_H
