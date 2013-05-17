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

#include <electronic/FluidSolver.h>
#include <electronic/Everything.h>
#include <electronic/NonlocalPCM.h>
#include <electronic/NonlinearPCM.h>
#include <electronic/LinearPCM.h>
#include <electronic/ConvCoupling.h>
#include <electronic/VDWCoupling.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/FluidMixture.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>
#include <fluid/IdealGasMonoatomic.h>
#include <core/DataIO.h>
#include <core/Units.h>


//--------------------------------------------------------------------------
//-------------------- Convolution coupled fluid solver --------------------
//--------------------------------------------------------------------------

//Extend FluidMixture to provide per-iteration dump option
class FluidMixtureJDFT : public FluidMixture
{
public:
	FluidMixtureJDFT(const Everything& e, double T)
	: FluidMixture(e.gInfo, T), e(e)
	{
	}
	
	bool report(int iter)
	{	//Call base class report:
		bool stateModified = FluidMixture::report(iter);
		//Dump, if required:
		((Dump&)e.dump)(DumpFreq_Fluid, iter);
		return stateModified;
	}
private:
	const Everything& e;
};


class ConvolutionJDFT : public FluidSolver
{
	FluidMixtureJDFT* fluidMixture;
	std::shared_ptr<ConvCoupling> coupling;
	std::shared_ptr<VDWCoupling> vdwCoupling;
	
	EnergyComponents Adiel; //fluid free energy components
	DataGptr Adiel_rhoExplicitTilde; //cached gradient of free energy of fluidMixture wrt rhoExplicit
	DataGptrCollection Ntilde; //cached densities of fluid

public:
	ConvolutionJDFT(const Everything& e, const FluidSolverParams& fsp)
	: FluidSolver(e, fsp), Adiel_rhoExplicitTilde(0)
	{
		//Initialize fluid mixture:
		fluidMixture = new FluidMixtureJDFT(e, fsp.T);
		fluidMixture->verboseLog = fsp.verboseLog;
		
		//Add the fluid components:
		for(const auto& c: fsp.components)
			c->addToFluidMixture(fluidMixture);

		fluidMixture->initialize(fsp.P, epsBulk, epsInf);

		//set fluid exCorr
		logPrintf("\n------- Fluid Exchange Correlation functional -------\n");
		((ExCorr&)fsp.exCorr).setup(e);
		
		//Initialize coupling:
		coupling = std::make_shared<ConvCoupling>(fluidMixture, fsp.exCorr);

		//Create van der Waals mixing functional
		assert(e.vanDerWaals);
		vdwCoupling = std::make_shared<VDWCoupling>(fluidMixture, e.vanDerWaals,
			e.vanDerWaals->getScaleFactor(fsp.exCorr.getName(), fsp.vdwScale));
	}

	~ConvolutionJDFT()
	{	delete fluidMixture;
	}

	bool needsGummel() { return true; }

	void loadState(const char* filename)
	{	fluidMixture->loadState(filename);
	}

	void saveState(const char* filename) const
	{	fluidMixture->saveState(filename);
	}

	void dumpDensities(const char* filenamePattern) const
	{	
		DataRptrCollection N; char filename[256];
		FluidMixture::Outputs outputs(&N,0,0);
		fluidMixture->getFreeEnergy(outputs);
		
		for(const auto& c: fsp.components)
			for(unsigned j=0; j<c->molecule.sites.size(); j++)
			{	const Molecule::Site& s = *(c->molecule.sites[j]);
				ostringstream oss; oss << "N_" << c->molecule.name;
				if(c->molecule.sites.size()>1) oss << "_" << s.name;
				sprintf(filename, filenamePattern, oss.str().c_str());
				logPrintf("Dumping %s... ", filename); logFlush();
				saveRawBinary(N[c->offsetDensity+j], filename);
				logPrintf("Done.\n"); logFlush();
			}
	}

	void dumpDebug(const char* filenamePattern) const
	{	
		DataRptrCollection N; char filename[256];
		FluidMixture::Outputs outputs(&N,0,0);
		fluidMixture->getFreeEnergy(outputs);
		
		//Save sphericalized site densities
		for(const auto& c: fsp.components)
			for(unsigned j=0; j<c->molecule.sites.size(); j++)
			{	const Molecule::Site& s = *(c->molecule.sites[j]);
				ostringstream oss; oss << "Nspherical_" << c->molecule.name;
				if(c->molecule.sites.size()>1) oss << "_" << s.name;
				sprintf(filename, filenamePattern, oss.str().c_str());
				logPrintf("Dumping %s... ", filename); logFlush();
				saveSphericalized(&N[c->offsetDensity+j], 1, filename);
				logPrintf("Done.\n"); logFlush();
			}
		
		//Adiel components
		string fname(filenamePattern);
		fname.replace(fname.find("%s"), 2, "Debug");
		logPrintf("Dumping '%s'... \t", fname.c_str());  logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());	
		fprintf(fp, "\nComponents of Adiel:\n");
		Adiel.print(fp, true, "   %18s = %25.16lf\n");	
		fclose(fp);
		
		//--------- if any, add additional explicit fluid debug output here!
	}

	void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
	{	
		//Set nCavity for nonlinear coupling functional
		coupling->setExplicit(nCavityTilde);
		//set rhoExplicit for electrostatic coupling
		fluidMixture->rhoExternal = clone(rhoExplicitTilde);
		if(!fluidMixture->state.size()) fluidMixture->initState(0.15, -3*fsp.T);
		if(!Adiel_rhoExplicitTilde) updateCached();
	}
	
	void updateCached()
	{	DataRptrCollection N;
		FluidMixture::Outputs outputs(&N, 0, &Adiel_rhoExplicitTilde, 0, &Adiel);
		
		fluidMixture->getFreeEnergy(outputs); //Fluid free energy including coupling
		Ntilde.resize(N.size());
		for(unsigned i=0; i<N.size(); i++)
			Ntilde[i] = J(N[i]);
	}

	void minimizeFluid()
	{	TIME("Fluid minimize", globalLog,
			fluidMixture->minimize(e.fluidMinParams);
			updateCached();
		)
	}
	
	double get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
	{
		assert(this->Adiel_rhoExplicitTilde); //Ensure that set() was called before calling get_Adiel_and_grad()
		Adiel_rhoExplicitTilde = clone(this->Adiel_rhoExplicitTilde);
		Adiel_nCavityTilde = 0; //clear previous, accumulate below
		extraForces.init(e.iInfo);
		
		//Update components of energy that depend on electronic state:
		EnergyComponents& Adiel = ((ConvolutionJDFT*)this)->Adiel;
		Adiel["ExtCoulomb"] = dot(fluidMixture->rhoExternal, O(Adiel_rhoExplicitTilde));
		Adiel["Fmix("+coupling->getName()+")"] = coupling->energyAndGrad(Ntilde, 0, &Adiel_nCavityTilde);
		Adiel["Fmix("+vdwCoupling->getName()+")"] = vdwCoupling->energyAndGrad(Ntilde, 0, &extraForces);
		return double(Adiel);
	}
};



//---------------------------------------------------------------------
//----------------  Interface to the electronic code ------------------
//---------------------------------------------------------------------

FluidSolver::FluidSolver(const Everything& e, const FluidSolverParams& fsp) : e(e), fsp(fsp)
{	//Initialize radial kernels in molecule sites:
	for(const auto& c: fsp.components)
		if(!c->molecule)
			c->molecule.setup(e.gInfo, c->Rvdw);

	//Set bulk dielectric constant
	if(fsp.epsBulkOverride)
		epsBulk = fsp.epsBulkOverride;
	else
	{	epsBulk = 1.;
		for(const auto& c: fsp.components)
			epsBulk += (c->Nbulk/c->pureNbulk(fsp.T)) * (c->epsBulk - 1.);
	}
	if(fsp.epsInfOverride)
		epsInf = fsp.epsInfOverride;
	else
	{	epsInf = 1.;
		for(const auto& c: fsp.components)
			epsInf += (c->Nbulk/c->pureNbulk(fsp.T)) * (c->epsInf - 1.);
	}

	//Check bulk charge balance and get screening prefactor:
	double NQ = 0., NQ2 = 0.;
	for(const auto& c: fsp.components)
	{	double Qmol = c->molecule.getCharge();
		if(Qmol && !c->Nnorm)
		{	NQ += c->Nbulk * Qmol;
			NQ2 += c->Nbulk * Qmol*Qmol;
		}
	}

	const double Qtol = 1e-12;
	if(fabs(NQ)>Qtol)
		die("Bulk fluid is non-neutral with a net charge density of %le e/bohr^3\n", NQ);
	k2factor = NQ2>Qtol ? (4*M_PI/fsp.T) * NQ2 : 0.;
}

FluidSolver* createFluidSolver(const Everything& e, const FluidSolverParams& fsp)
{	if(fsp.fluidType != FluidNone)
		Citations::add("Framework of Joint Density Functional Theory", "S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)");
	logPrintf("%s", fsp.initWarnings.c_str());
	switch(fsp.fluidType)
	{	case FluidNone:
			return 0; //No solver needed
		case FluidLinearPCM:
			return new LinearPCM(e, fsp);
		case FluidNonlinearPCM:
			return new NonlinearPCM(e, fsp);
		case FluidNonlocalPCM:
			return new NonlocalPCM(e, fsp);
		case FluidClassicalDFT:
			return new ConvolutionJDFT(e, fsp);
		default:
			assert(!"Unknown fluid type");
			return 0;
	}
}
