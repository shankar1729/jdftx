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

#include <fluid/FluidSolver.h>
#include <electronic/Everything.h>
#include <fluid/SaLSA.h>
#include <fluid/NonlinearPCM.h>
#include <fluid/LinearPCM.h>
#include <fluid/ConvCoupling.h>
#include <fluid/VDWCoupling.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/FluidMixture.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>
#include <fluid/Fex_LJ.h>
#include <fluid/IdealGasMonoatomic.h>
#include <core/ScalarFieldIO.h>
#include <core/Units.h>
#include <core/WignerSeitz.h>

//--------------------------------------------------------------------------
//-------------------- Convolution coupled fluid solver --------------------
//--------------------------------------------------------------------------

//Extend FluidMixture to provide per-iteration dump option
class FluidMixtureJDFT : public FluidMixture
{
public:
	FluidMixtureJDFT(const Everything& e, const GridInfo& gInfo, double T)
	: FluidMixture(gInfo, T), e(e)
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
	std::vector<std::shared_ptr<Fmix>> FmixPtr;
	
	EnergyComponents Adiel; //fluid free energy components
	ScalarFieldTilde Adiel_rhoExplicitTilde; //cached gradient of free energy of fluidMixture wrt rhoExplicit
	ScalarFieldTildeArray Ntilde; //cached densities of fluid

public:
	ConvolutionJDFT(const Everything& e, const FluidSolverParams& fsp)
	: FluidSolver(e, fsp), Adiel_rhoExplicitTilde(0)
	{
		//Initialize fluid mixture:
		fluidMixture = new FluidMixtureJDFT(e, gInfo, fsp.T);
		fluidMixture->verboseLog = fsp.verboseLog;
		
		//Add the fluid components:
		for(const auto& c: fsp.components)
			c->addToFluidMixture(fluidMixture);

		if(fsp.FmixList.size())
		{
		        //create fluid mixtures
		        logPrintf("\n------------ Fluid Mixing Functionals ------------\n");
			for(const auto& f: fsp.FmixList)
			{
			       std::shared_ptr<FluidComponent> c1 = f.fluid1;
			       string name1 = c1->molecule.name;
			       std::shared_ptr<FluidComponent> c2 = f.fluid2;
			       string name2 = c2->molecule.name;
		      
			       std::shared_ptr<Fmix> Fmix;
			       if (f.FmixType == GaussianKernel)
				 Fmix = std::make_shared<Fmix_GaussianKernel>(fluidMixture,c1,c2,f.energyScale,f.lengthScale);	
			       else if (f.FmixType == LJPotential)
				 Fmix = std::make_shared<Fmix_LJ>(fluidMixture,c1,c2,f.energyScale,f.lengthScale);
			       else
				 die("Valid mixing functional between %s and %s not specified!\n",name1.c_str(),name2.c_str());
			       FmixPtr.push_back(Fmix);		      
			}
		}

		fluidMixture->initialize(fsp.P, epsBulk, epsInf);

		//set fluid exCorr
		logPrintf("\n------- Fluid Exchange Correlation functional -------\n");
		((ExCorr&)fsp.exCorr).setup(e);
		
		//Initialize coupling:
		coupling = std::make_shared<ConvCoupling>(fluidMixture, fsp.exCorr);

		//Create van der Waals mixing functional
		assert(e.vanDerWaals);
		vdwCoupling = std::make_shared<VDWCoupling>(fluidMixture, atpos, e.vanDerWaals,
			e.vanDerWaals->getScaleFactor(fsp.exCorr.getName(), fsp.vdwScale));

		//---- G=0 constraints -----

		//Electron density in the bulk of the fluid
		double nFl_bulk = 0.0;	

		for(const auto& c: fsp.components) 
		{  
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
			{
			  const Molecule::Site& s = *(c->molecule.sites[i]);
			  nFl_bulk += c->idealGas->get_Nbulk()*s.elecKernel(0)*s.positions.size();
			}
		}
		
		logPrintf("\nBulk electron density of the liquid: %le bohr^-3\n",nFl_bulk);

		//calculate G=0 offset due to coupling functional evaluated at bulk fluid density
		ScalarField nBulk;
		nullToZero(nBulk,e.gInfo);
		//initialize constant ScalarField with density nFl_bulk
		for (int i=0; i<e.gInfo.nr; i++)
			nBulk->data()[i] = nFl_bulk;
		ScalarField Vxc_bulk;
		(coupling->exCorr)(nBulk, &Vxc_bulk, true);
		logPrintf("Electron deep in fluid experiences coupling potential: %lg H\n\n", Vxc_bulk->data()[0]);
		coupling->Vxc_bulk = Vxc_bulk->data()[0];
	}

	~ConvolutionJDFT()
	{	delete fluidMixture;
	}

	bool prefersGummel() const { return true; }

	double bulkPotential() {return coupling->Vxc_bulk;}

	void loadState(const char* filename)
	{	fluidMixture->loadState(filename);
	}

	void saveState(const char* filename) const
	{	fluidMixture->saveState(filename);
	}

	void dumpDensities(const char* filenamePattern) const
	{	
		ScalarFieldArray N; char filename[256];
		FluidMixture::Outputs outputs(&N,0,0);
		fluidMixture->getFreeEnergy(outputs);
		
		for(const auto& c: fsp.components)
			for(unsigned j=0; j<c->molecule.sites.size(); j++)
			{	const Molecule::Site& s = *(c->molecule.sites[j]);
				ostringstream oss; oss << "N_" << c->molecule.name;
				if(c->molecule.sites.size()>1) oss << "_" << s.name;
				sprintf(filename, filenamePattern, oss.str().c_str());
				logPrintf("Dumping %s... ", filename); logFlush();
				if(mpiWorld->isHead()) saveRawBinary(N[c->offsetDensity+j], filename);
				logPrintf("Done.\n"); logFlush();
			}
	}

	void dumpDebug(const char* filenamePattern) const
	{	
		ScalarFieldArray N; char filename[256];
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
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());	
			fprintf(fp, "\nComponents of Adiel:\n");
			Adiel.print(fp, true, "   %18s = %25.16lf\n");	
			fclose(fp);
		}
		logPrintf("Done.\n"); logFlush();
		
		//--------- if any, add additional explicit fluid debug output here!
	}

	void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
	{	
		//Set nCavity for nonlinear coupling functional
		coupling->setExplicit(nCavityTilde);
		//set rhoExplicit for electrostatic coupling
		fluidMixture->rhoExternal = rhoExplicitTilde;
		if(!fluidMixture->state.size()) fluidMixture->initState(0.15, -3*fsp.T);
		if(!Adiel_rhoExplicitTilde) updateCached();
	}
	
	void updateCached()
	{	ScalarFieldArray N;
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
	
	double get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces) const
	{
		assert(this->Adiel_rhoExplicitTilde); //Ensure that set() was called before calling get_Adiel_and_grad()
		Adiel_rhoExplicitTilde = clone(this->Adiel_rhoExplicitTilde);
		Adiel_nCavityTilde = 0; //clear previous, accumulate below
		if(extraForces) extraForces->init(e.iInfo);
		
		//Update components of energy that depend on electronic state:
		EnergyComponents& Adiel = ((ConvolutionJDFT*)this)->Adiel;
		Adiel["ExtCoulomb"] = dot(fluidMixture->rhoExternal, O(Adiel_rhoExplicitTilde));
		Adiel["Fmix("+coupling->getName()+")"] = coupling->energyAndGrad(Ntilde, 0, &Adiel_nCavityTilde);
		Adiel["Fmix("+vdwCoupling->getName()+")"] = vdwCoupling->energyAndGrad(Ntilde, 0, extraForces);
		return double(Adiel);
	}
};



//---------------------------------------------------------------------
//----------------  Interface to the electronic code ------------------
//---------------------------------------------------------------------

FluidSolver::FluidSolver(const Everything& e, const FluidSolverParams& fsp)
: e(e), gInfo(e.coulomb->gInfo), fsp(fsp), atpos(e.iInfo.species.size())
{	//Initialize radial kernels in molecule sites:
	for(const auto& c: fsp.components)
		if(!c->molecule)
			c->molecule.setup(gInfo, c->Rvdw);

	//Set bulk dielectric constant
	if(fsp.epsBulkTensor.length_squared()) //set epsBulk to be consistent with anisotropic override:
		epsBulk = (1./3)*(fsp.epsBulkTensor[0] + fsp.epsBulkTensor[1] + fsp.epsBulkTensor[2]);
	else if(fsp.epsBulkOverride)
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
	if(fsp.screenOverride)
		k2factor = epsBulk / pow(fsp.screenOverride, 2);
	else
		k2factor = NQ2>Qtol ? (4*M_PI/fsp.T) * NQ2 : 0.;
	
	if(e.iInfo.ionWidth)
		logPrintf("\nCorrection to mu due to finite nuclear width = %lg\n", ionWidthMuCorrection());
}

double FluidSolver::ionWidthMuCorrection() const
{	return (4*M_PI/gInfo.detR) * (-0.5*pow(e.iInfo.ionWidth,2)) * e.iInfo.getZtot();
}

void FluidSolver::set(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{	for(unsigned iSp=0; iSp<atpos.size(); iSp++)
		atpos[iSp] = e.iInfo.species[iSp]->atpos;
	if(e.coulombParams.embed)
	{	matrix3<> embedScaleMat = Diag(e.coulomb->embedScale); //lattice coordinate scale factor due to embedding
		for(std::vector<vector3<> >& posArr: atpos)
			for(vector3<>& pos: posArr) //transform to embedded lattice coordinates:
				pos = embedScaleMat *  e.coulomb->wsOrig->restrict(pos - e.coulomb->xCenter);
		ScalarFieldTilde rhoExplicitTildeExpand = e.coulomb->embedExpand(rhoExplicitTilde);
		if(!k2factor) rhoExplicitTildeExpand->setGzero(0.); //No screening => apply neutralizing background charge
		set_internal(rhoExplicitTildeExpand, e.coulomb->embedExpand(nCavityTilde));
	}
	else
	{	if(!k2factor) ((ScalarFieldTilde&)rhoExplicitTilde)->setGzero(0.); //No screening => apply neutralizing background charge
		set_internal(rhoExplicitTilde, nCavityTilde);
	}
}

double FluidSolver::get_Adiel_and_grad(ScalarFieldTilde* Adiel_rhoExplicitTilde, ScalarFieldTilde* Adiel_nCavityTilde, IonicGradient* extraForces) const
{	if(e.coulombParams.embed)
	{	ScalarFieldTilde Adiel_rho_big, Adiel_n_big;
		double Adiel = get_Adiel_and_grad_internal(Adiel_rho_big, Adiel_n_big, extraForces);
		if(A_rhoNonES) ((FluidSolver*)this)->A_rhoNonES = e.coulomb->embedShrink(A_rhoNonES);
		if(Adiel_rhoExplicitTilde) *Adiel_rhoExplicitTilde = e.coulomb->embedShrink(Adiel_rho_big);
		if(Adiel_nCavityTilde) *Adiel_nCavityTilde = e.coulomb->embedShrink(Adiel_n_big);
		if(extraForces) *extraForces = Diag(e.coulomb->embedScale) * (*extraForces); //transform to original contravariant lattice coordinates
		return Adiel;
	}
	else
	{	ScalarFieldTilde Adiel_rho_temp, Adiel_n_temp;
		return get_Adiel_and_grad_internal(
			Adiel_rhoExplicitTilde ? (*Adiel_rhoExplicitTilde) : Adiel_rho_temp,
			Adiel_nCavityTilde ? (*Adiel_nCavityTilde) : Adiel_n_temp, extraForces);
	}
}

void FluidSolver::getSusceptibility(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldTildeArray& sTilde, bool elecOnly) const
{	ScalarFieldArray sArr;
	getSusceptibility_internal(omega, susceptibility, sArr, elecOnly);
	sTilde.clear();
	for(const ScalarField& s: sArr)
	{	if(e.coulombParams.embed)
			sTilde.push_back(e.coulomb->embedShrink(J(s)));
		else
			sTilde.push_back(J(s));
	}
}

void FluidSolver::getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const
{	die("\nSusceptibility not yet implemented for this fluid type.\n\n");
}


FluidSolver* createFluidSolver(const Everything& e, const FluidSolverParams& fsp)
{	logPrintf("%s", fsp.initWarnings.c_str());
	switch(fsp.fluidType)
	{	case FluidNone:
			return 0; //No solver needed
		case FluidLinearPCM:
			return new LinearPCM(e, fsp);
		case FluidNonlinearPCM:
			return new NonlinearPCM(e, fsp);
		case FluidSaLSA:
			return new SaLSA(e, fsp);
		case FluidClassicalDFT:
			return new ConvolutionJDFT(e, fsp);
		default:
			assert(!"Unknown fluid type");
			return 0;
	}
}
