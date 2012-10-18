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
#include <fluid/IdealGasMuEps.h>
#include <fluid/FluidMixture.h>
#include <fluid/Fex_H2O_Lischner10.h>
#include <fluid/Fex_H2O_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>
#include <fluid/Fex_HardSphereIon.h>
#include <fluid/IdealGasMonoatomic.h>
#include <fluid/Fmix_IonSolvation.h>
#include <core/DataIO.h>
#include <core/Units.h>


//--------------------------------------------------------------------------
//-------------------- Convolution coupled fluid solver --------------------
//--------------------------------------------------------------------------

//Constrain center of cell to be at 0 electrostatic potential
double zeroCenter(const DataRptr& dtot, DataRptr& dd0_ddtot)
{
	const GridInfo& gInfo = dtot->gInfo;
	initZero(dd0_ddtot, gInfo);
	int centerIndex = gInfo.fullRindex(vector3<int>(gInfo.S[0]/2, gInfo.S[1]/2, gInfo.S[2]/2));
	double d0 = -dtot->data()[centerIndex];
	dd0_ddtot->data()[centerIndex] = -1.0;
	return d0;
}


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
		if(e.dump.shouldDump(DumpFreq_Fluid, iter))
			((Dump&)e.dump)(DumpFreq_Fluid);
		return stateModified;
	}
private:
	const Everything& e;
};


class ConvolutionJDFT : public FluidSolver
{
	FluidMixtureJDFT* fluidMixture;
	Fex* fex;
	
	SO3quad* quad;
	TranslationOperatorSpline* trans;
	IdealGasMuEps* idgas;
	
	ConvCoupling* coupling;
	
	double* AwaterCached; //cached free energy of fluidMixture (without the coupling part)
	DataGptr* grad_rhoExplicitTildeCached; //cached gradient of free energy of fluidMixture wrt rhoExplicit
	
public:
	ConvolutionJDFT(const Everything& e, FluidSolverParams& params, FluidType type)
	: FluidSolver(e), AwaterCached(0), grad_rhoExplicitTildeCached(0)
	{
		
		//Initialize fluid mixture:
		fluidMixture = new FluidMixtureJDFT(e, params.T);
		fluidMixture->verboseLog = params.verboseLog;
		
		//Initialize excess functional:
		switch(type)
		{	
			case FluidLischner10: fex = new Fex_H2O_Lischner10(*fluidMixture); break;
			case FluidScalarEOS: fex = new Fex_H2O_ScalarEOS(*fluidMixture); break;
			case FluidBondedVoids: fex = new Fex_H2O_BondedVoids(*fluidMixture); break;
			case FluidHSIonic: break;
			default: assert(!"This is not a JDFT3 functional");
		}
		
		if (type != FluidHSIonic)
		{
			//Initialize ideal gas:
			const int Zn = 2; //Water molecule has Z2 symmetry about dipole axis
			quad = new SO3quad(params.s2quadType, Zn, params.quad_nBeta, params.quad_nAlpha, params.quad_nGamma);
			trans = new TranslationOperatorSpline(e.gInfo, TranslationOperatorSpline::Linear);
			idgas = new IdealGasMuEps(fex, 55.0*mol/liter, *quad, *trans);
		}
		
		int nIons = params.hSIons.size();

		if((type == FluidHSIonic) && nIons == 0)	
			die("At least one fluid-ion must be specified in FluidHSIonic.");
		
		for(int iIon = 0; iIon < nIons; iIon++)
		{                                 
			HardSphereIon* Ion = &params.hSIons[iIon];
			Ion->fex = new Fex_HardSphereIon(*fluidMixture, Ion);
			Ion->idgas = new IdealGasMonoatomic(Ion->fex, Ion->Concentration);
			if (Ion->MixFunc)
				Ion->fmix = new Fmix_IonSolvation(*fluidMixture, *fex, Ion);                                                                   
		}	
		
		//Initialize coupling:
		coupling = new ConvCoupling(*fluidMixture);

		//set fluid exCorr
		logPrintf("\n------- Fluid Exchange Correlation functional -------\n");
		params.exCorr.setup(e);
		coupling->exCorr = &params.exCorr;
		
		
		//Set the electronic site density model for the coupling -- specific for scalarEOS water
		FluidMixture::Component& water = (FluidMixture::Component&) fluidMixture->get_component(0);
		
		//Set the electronic site density model for the coupling
		//For now set to be the bad model used as default in the older versions of JDFTx 	
		SiteProperties& Osite = *water.indexedSite[0];
		Osite.siteName="O";
		Osite.couplingZnuc = params.oxygenZnuc;
						
		SiteProperties& Hsite = *water.indexedSite[1];
		Hsite.siteName="H";
		Hsite.couplingZnuc = params.hydrogenZnuc;
		
		switch(params.convCouplingH2OModel)
		{	
			case ConvCouplingExponential: 
			{		
				Osite.convCouplingSiteCharge = params.oxygenSiteCharge;
				Osite.convCouplingWidth = params.oxygenWidth;
				Osite.kernelFilename = params.oxygenFilename;
				coupling->setExponentialKernel(Osite);
				
				if(Osite.kernelFilename.length()!=0)
					coupling->setRadialKernel(Osite);
				
				Hsite.convCouplingSiteCharge = params.hydrogenSiteCharge;
				if(params.hydrogenSiteCharge+Osite.convCouplingSiteCharge/2.0>1e-12)
					die("Water molecule has net charge due to unbalanced site charges in fluid coupling.\n");
				Hsite.convCouplingWidth = params.hydrogenWidth;
				Hsite.kernelFilename = params.hydrogenFilename;
				coupling->setExponentialKernel(Hsite);
				
				if(Hsite.kernelFilename.length()!=0)
					coupling->setRadialKernel(Hsite);
				
				break;
				
			}
			case ConvCouplingExpCuspless: 
			{		
				Osite.convCouplingSiteCharge = params.oxygenSiteCharge;
				Osite.convCouplingWidth = params.oxygenWidth;
				Osite.kernelFilename = params.oxygenFilename;
				coupling->setExpCusplessKernel(Osite);
				
				if(Osite.kernelFilename.length()!=0)
					coupling->setRadialKernel(Osite);
				
				Hsite.convCouplingSiteCharge = params.hydrogenSiteCharge;
				if(params.hydrogenSiteCharge+Osite.convCouplingSiteCharge/2.0>1e-12)
					die("Water molecule has net charge due to unbalanced site charges in fluid coupling.\n");
				Hsite.convCouplingWidth = params.hydrogenWidth;
				Hsite.kernelFilename = params.hydrogenFilename;
				coupling->setExpCusplessKernel(Hsite);
				
				if(Hsite.kernelFilename.length()!=0)
					coupling->setRadialKernel(Hsite);
				
				break;
				
			}
			case ConvCouplingRadialFunction:
			{	
				Osite.kernelFilename = params.oxygenFilename;
				Hsite.kernelFilename = params.hydrogenFilename;
				
				coupling->setRadialKernel(Osite);
				coupling->setRadialKernel(Hsite);
				
				break;
			}
			case ConvCouplingBinaryKernel:
			{	
				Osite.kernelFilename = params.oxygenFilename;
				Hsite.kernelFilename = params.hydrogenFilename;
				
				coupling->setBinaryKernel(Osite);
				coupling->setBinaryKernel(Hsite);
				break;
			}
			case ConvCouplingNone: break;
			default:
			{
				die("Unknown convolution coupling model specified for water.\n")
			}
			
		}
		
		for(int iIon = 0; iIon < nIons; iIon++)
		{  
			HardSphereIon* Ion = &params.hSIons[iIon];
			FluidMixture::Component& IonComp = (FluidMixture::Component&) fluidMixture->get_component(iIon+1);
			if (Ion->convCouplingModel!=ConvCouplingNone)
			{
				SiteProperties& IonSite = *IonComp.indexedSite[0];
				IonSite.siteName = Ion->name;
				IonSite.convCouplingSiteCharge = Ion->Z;
				
				switch(Ion->convCouplingModel)
				{
					case ConvCouplingExponential:
					{
						IonSite.couplingZnuc = Ion->Znuc;
						IonSite.convCouplingWidth = Ion->CouplingWidth;
						IonSite.kernelFilename = Ion->CouplingFilename;
						coupling->setExponentialKernel(IonSite);
						
						if(IonSite.kernelFilename.length()!=0)
						coupling->setRadialKernel(IonSite);
						break;
					}
					case ConvCouplingExpCuspless:
					{
						IonSite.couplingZnuc = Ion->Znuc;
						IonSite.convCouplingWidth = Ion->CouplingWidth;
						IonSite.kernelFilename = Ion->CouplingFilename;
						coupling->setExpCusplessKernel(IonSite);
						
						if(IonSite.kernelFilename.length()!=0)
						coupling->setRadialKernel(IonSite);
						break;
					}
					case ConvCouplingRadialFunction:
					{
						IonSite.kernelFilename = Ion->CouplingFilename;
						coupling->setRadialKernel(IonSite);
						break;
					}
					case ConvCouplingBinaryKernel:
					{ 
						IonSite.kernelFilename = Ion->CouplingFilename;
						coupling->setBinaryKernel(IonSite);
						break;
					}	
					case ConvCouplingNone: break;
					default:
					{
						die("Unknown convolution coupling model specified for ion %s.\n",IonSite.siteName.c_str())
					}
				}
			}
		}
		
		fluidMixture->setPressure(1.01325*Bar);
		
		//---- G=0 constraint -----
		fluidMixture->d0calc = zeroCenter;
	}

	~ConvolutionJDFT()
	{	
		delete coupling;
		delete idgas;
		delete trans;
		delete quad;
		delete fex;
		delete fluidMixture;
		
		if(AwaterCached) delete AwaterCached;
		if(grad_rhoExplicitTildeCached) delete grad_rhoExplicitTildeCached;
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
		
		//loop over sites and use the chain rule to calculate nonlinear contribution to the coupling.
		for(unsigned ic=0; ic<fluidMixture->get_nComponents(); ic++)
		{		
			const FluidMixture::Component& c = fluidMixture->get_component(ic);
			string moleculeName = c.molecule->name;
			if (c.molecule->nIndices == 1)
			{
				sprintf(filename, "%sN_%s", filenamePattern, moleculeName.c_str());
				logPrintf("Dumping %s... ", filename); logFlush();
				saveRawBinary(N[c.offsetDensity], filename); logPrintf("Done.\n"); logFlush();
			}
			else
			{
				for(int j=0; j<c.molecule->nIndices; j++)
				{
					const SiteProperties& s = *c.indexedSite[j];
					sprintf(filename, "%sN_%s_%s", filenamePattern, moleculeName.c_str(), s.siteName.c_str());
					logPrintf("Dumping %s... ", filename); logFlush();
					saveRawBinary(N[c.offsetDensity+j], filename); logPrintf("Done.\n"); logFlush();
				}
			}
		}	
	}

	void dumpDebug(const char* filenamePattern) const
	{	/*
		coupling.dumpDebug(filenamePattern);
		//Compute potential on O and H sites (done here instead of inside coupling, since need access to densities)
		DataRptrOH N, grad_N; fluidMixture->getFreeEnergy(&N); EnergyComponents omega;
		DataGptrOH Ntilde=J(N), grad_Ntilde; DataRptrVec p, grad_p, eps, grad_eps;
		MomentFields mf = {omega, N, Ntilde, p, eps, grad_N, grad_Ntilde, grad_p, grad_eps};
		nullToZero(grad_N, e.gInfo); nullToZero(grad_Ntilde, e.gInfo);
		((ConvolutionJDFT*)this)->coupling.calculateCoupling(mf); //and now grad_N and grad_Ntilde contain the coupling potential!
		grad_N += Jdag(grad_Ntilde); grad_Ntilde={0,0}; //now all the coupling potential is in real space

		char filename[256];
		sprintf(filename, filenamePattern, "VO");
		logPrintf("Dumping %s... ", filename); logFlush();
		saveRawBinary(grad_N.O(), filename); logPrintf("Done.\n"); logFlush();
		sprintf(filename, filenamePattern, "VH");
		logPrintf("Dumping %s... ", filename); logFlush();
		saveRawBinary(grad_N.H(), filename); logPrintf("Done.\n"); logFlush();
		*/
		//--------- if any, add additional JDFT3 debug output here!
	}

	void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
	{	
		//Set nCavity for nonlinear coupling functional
		coupling->setExplicit(nCavityTilde);
		//set rhoExplicit for electrostatic coupling
		fluidMixture->rhoExternal=clone(rhoExplicitTilde);
		if((fluidMixture->state.size())==0) fluidMixture->initState(0.05);
		if(!AwaterCached)
		{	
			AwaterCached = new double;
			*AwaterCached=0.0;
			grad_rhoExplicitTildeCached = new DataGptr;
			nullToZero(*grad_rhoExplicitTildeCached,e.gInfo);
			
			FluidMixture::Outputs outputs(0,0,grad_rhoExplicitTildeCached);
			
			*AwaterCached = fluidMixture->getFreeEnergy(outputs); //Fluid free energy including coupling
			*AwaterCached -= dot(rhoExplicitTilde,O(*grad_rhoExplicitTildeCached)); //subtract Electrostatic coupling energy
			*AwaterCached -= coupling->computeElectronic(); //subtract Nonlinear coupling energy
			
		}
		
	}

	void minimizeFluid()
	{	
		assert(AwaterCached); //Ensure that set() was called before calling minimize_fluid()
		TIME("Fluid minimize", globalLog,
			fluidMixture->minimize(e.fluidMinParams);
			
			FluidMixture::Outputs outputs(0,0,grad_rhoExplicitTildeCached);
			
			*AwaterCached = fluidMixture->getFreeEnergy(outputs); //Fluid free energy including coupling
			*AwaterCached -= dot(fluidMixture->rhoExternal,O(*grad_rhoExplicitTildeCached)); //subtract Electrostatic coupling energy
			*AwaterCached -= coupling->computeElectronic(); //subtract Nonlinear coupling energy
		)
	}

	double get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde)
	{
		assert(AwaterCached); //Ensure that set() was called before calling get_Adiel_and_grad()
		
		grad_rhoExplicitTilde = clone(*grad_rhoExplicitTildeCached);	
		//Uncomment to print the components of Adiel (fluid alone, electrostatic coupling, convolution coupling) 
		//logPrintf("Awater: %le ExtCoulomb: %le Coupling: %le\n", *AwaterCached, 
		//		dot(fluidMixture->rhoExternal,O(*grad_rhoExplicitTildeCached)), coupling->computeElectronic());
		return *AwaterCached + dot(fluidMixture->rhoExternal,O(*grad_rhoExplicitTildeCached))
				+ coupling->computeElectronic(&grad_nCavityTilde);
	}
};



//---------------------------------------------------------------------
//----------------  Interface to the electronic code ------------------
//---------------------------------------------------------------------

FluidSolver::FluidSolver(const Everything& e) : e(e)
{
}

FluidSolver* createFluidSolver(FluidType type, const Everything& e, FluidSolverParams& params)
{	switch(type)
	{
		case FluidNone: return 0; //No solver needed
		case FluidLinear: return new LinearPCM(e, params);
		case FluidLinearPCM: return new LinearPCM(e, params);
		case FluidNonlinearPCM: return new NonlinearPCM(e, params);
		case FluidNonlocalPCM: return new NonlocalPCM(e, params);
		default: //All JDFT3 functionals:
			return new ConvolutionJDFT(e, params, type);
	}
}
