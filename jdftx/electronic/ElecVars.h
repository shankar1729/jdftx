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

#include <fluid/FluidSolverParams.h>
#include <core/ScalarFieldArray.h>
#include <string>
#include <memory>

struct ElecGradient;

//! @addtogroup ElecSystem
//! @{

//! Electronic variables and main energy calculator
class ElecVars
{
public:
	std::vector<ColumnBundle> C; //!< orthonormal electronic wavefunctions
	std::vector<diagMatrix> Haux_eigs; //!< auxilliary hamiltonian eigenvalues
	std::vector<diagMatrix> F;  //!< the fillings (diagonal matrices) for each state
	
	std::vector<matrix> Hsub; //!< Subspace Hamiltonian:  Hsub[q]=C[q]^H*C[q]
	std::vector<matrix> Hsub_evecs; //!< eigenvectors of Hsub[q] in columns
	std::vector<diagMatrix> Hsub_eigs; //!< eigenvalues of Hsub[q]
	
	std::vector< std::vector<matrix> > VdagC; //!< cached pseudopotential projections (by state and then species)
	
	//Densities and potentials:
	ScalarFieldArray n; //!< electron density (single ScalarField) or spin density (two ScalarFields [up,dn]) or spin density matrix (four ScalarFields [UpUp, DnDn, Re(UpDn), Im(UpDn)])
	ScalarFieldArray nAccumulated; //!< ElecVars::n accumulated over an MD trajectory
	ScalarFieldArray get_nXC() const; //!< return the total (spin) density including core contributions
	ScalarField get_nTot() const { return n.size()==1 ? n[0] : n[0]+n[1]; } //!< return the total electron density (even in spin polarized situations)
	
	ScalarFieldArray tau; //!< kinetic energy density including tauCore, if present (computed if a meta GGA is being used)
	
	ScalarFieldTilde d_fluid; //!< electrostatic potential due to fluid
	ScalarFieldTilde V_cavity; //!< non-electrostatic potential on electrons due to fluid
	
	ScalarFieldArray Vscloc; //! Local part of (optionally spin-dependent) self-consistent potential
	ScalarFieldArray Vtau; //! Gradient w.r.t kinetic energy density (if meta-GGA)
	
	std::vector<matrix> rhoAtom, U_rhoAtom; //!< Atomic density matrices and gradients w.r.t them (for DFT+U)
	
	//External interactions:
	ScalarFieldArray Vexternal; //!< external potential
	ScalarFieldTilde rhoExternal; //!< external charge density
	bool rhoExternalSelfEnergy; //!< whether to include self-energy of rhoExternal in output energy
	struct BoxPotential //!< box potetial desciptor
	{	vector3<> min; //!< minimum of bounding box
		vector3<> max; //!< maximum of bounding box
		double Vin; //!< potential inside
		double Vout; //!< potential outside
		double convolve_radius; //!< smoothing radius
	};
	std::vector<BoxPotential> boxPot; //!< parameters for the external box potential

	//Fluid properties
	FluidSolverParams fluidParams;
	std::shared_ptr<struct FluidSolver> fluidSolver;
	string fluidInitialStateFilename;

	//Wavefunction initialization:
	string wfnsFilename; //!< file to read wavefunctions from
	std::shared_ptr<struct ColumnBundleReadConversion> readConversion; //!< ColumnBundle conversion
	bool isRandom; //!< indicates whether the electronic state is random (not yet minimized)
	bool initLCAO; //!< initialize wave functions using linear combinations of atomic orbitals
	bool skipWfnsInit; //!< whether to skip wavefunction initialization (used to speed up dry runs, phonon calculations)

	string eigsFilename; //!< file to read eigenvalues from
	
	//Auxiliary hamiltonian initialization
	bool HauxInitialized; //!< whether Haux has been read in/computed
	
	string nFilenamePattern; //!< file pattern to read electron (spin,kinetic) density from
	string VFilenamePattern; //!< file pattern to read electron (spin,kinetic) potential from
	
	ElecVars();
	void setup(const Everything &everything);

	//! Compute the terms written as a functional of the electronic density, and its gradient i.e. Vscloc
	//! If supplied, alternateExCorr replaces the main exchange and correlaton functional
	void EdensityAndVscloc(Energies& ener, const ExCorr* alternateExCorr=0);
	
	//! Update and return the electronic system energy.
	//! Optionally compute the gradient, preconditioned gradient and/or the subspace hamiltonian
	double elecEnergyAndGrad(Energies& ener, ElecGradient* grad=0, ElecGradient* Kgrad=0, bool calc_Hsub=false); 
	
	//! Set C to eigenvectors of the subspace hamiltonian
	void setEigenvectors(); 
	
	//! Compute the kinetic energy density
	ScalarFieldArray KEdensity() const;
	
	//! Calculate density using current orthonormal wavefunctions (C)
	ScalarFieldArray calcDensity() const;
	
	//! Orthonormalise wavefunctions, with an optional extra rotation
	//! If extraRotation is present, it is applied after symmetric orthononormalization,
	//! and on output extraRotation contains the net transformation applied to the wavefunctions.
	void orthonormalize(int q, matrix* extraRotation=0);
	
	//! Applies the Kohn-Sham Hamiltonian on the orthonormal wavefunctions C, and computes Hsub if necessary, for a single quantum number
	//! Returns the Kinetic energy contribution from q, which can be used for the inverse kinetic preconditioner
	double applyHamiltonian(int q, const diagMatrix& Fq, ColumnBundle& HCq, Energies& ener, bool need_Hsub = false);
	
private:
	const Everything* e;
	
	std::vector<string> VexternalFilename; //!< external potential filename (read in real space)
	friend struct CommandVexternal;
	
	string rhoExternalFilename; //!< external charge filename
	friend struct CommandRhoExternal;

	int lcaoIter; //!< number of iterations for LCAO (automatic if negative)
	double lcaoTol; //!< tolerance for LCAO subspace minimization
	int LCAO(); //!< Initialize LCAO wavefunctions (returns the number of bands initialized)
	friend struct CommandWavefunction;
	friend struct CommandLcaoParams;
	friend class Dump;
};

//! @}
#endif // JDFTX_ELECTRONIC_ELECVARS_H
