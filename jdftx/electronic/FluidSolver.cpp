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
		((Dump&)e.dump)(DumpFreq_Fluid, iter);
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
	VDWCoupling* vdwCoupling;
	
	double* AwaterCached; //cached free energy of fluidMixture (without the coupling part)
	DataGptr* Adiel_rhoExplicitTildeCached; //cached gradient of free energy of fluidMixture wrt rhoExplicit
	DataRptrCollection* NCached; //cached densities of fluid for use in van der Waals correction
	
public:
	ConvolutionJDFT(const Everything& e, FluidSolverParams& params)
	: FluidSolver(e), AwaterCached(0), Adiel_rhoExplicitTildeCached(0), NCached(0)
	{
		
		//Initialize fluid mixture:
		fluidMixture = new FluidMixtureJDFT(e, params.T);
		fluidMixture->verboseLog = params.verboseLog;
		
		//Initialize excess functional:
		switch(params.fluidType)
		{	case FluidFittedCorrelations: fex = new Fex_H2O_FittedCorrelations(*fluidMixture); break;
			case FluidScalarEOS: fex = new Fex_H2O_ScalarEOS(*fluidMixture); break;
			case FluidScalarEOSCustom: fex = new Fex_H2O_Custom(*fluidMixture, params.H2OSites); break;
			case FluidBondedVoids: fex = new Fex_H2O_BondedVoids(*fluidMixture); break;
			case FluidHSIonic: break;
			default: assert(!"This is not a JDFT3 functional");
		}
		
		
		if (params.fluidType != FluidHSIonic)
		{
			//Initialize ideal gas:
			const int Zn = 2; //Water molecule has Z2 symmetry about dipole axis
			quad = new SO3quad(params.s2quadType, Zn, params.quad_nBeta, params.quad_nAlpha, params.quad_nGamma);
			trans = new TranslationOperatorSpline(e.gInfo, TranslationOperatorSpline::Linear);
			idgas = new IdealGasMuEps(fex, 55.0*mol/liter, *quad, *trans);
		}
		
		int nIons = params.hSIons.size();

		if((params.fluidType == FluidHSIonic) && nIons == 0)	
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
			
		if (params.fluidType == FluidScalarEOSCustom) //for right now keep this separate to avoid creating bugs.
		{
			for (int iSite=0; iSite < water.molecule->nIndices; iSite++)
			{
				//Set the electronic site density model for the coupling
				SiteProperties& water_site = *water.indexedSite[iSite];
				
		
				switch(water_site.convCouplingSiteModel)
				{	
					case ConvCouplingExponential: 
					{		
						coupling->setExponentialKernel(water_site);
				
						if(water_site.kernelFilename.length()!=0)
							coupling->setRadialKernel(water_site);
				
						break;
					}
					case ConvCouplingExpCuspless: 
					{	
						coupling->setExpCusplessKernel(water_site);
				
						if(water_site.kernelFilename.length()!=0)
							coupling->setRadialKernel(water_site);
				
						break;				
					}
					case ConvCouplingRadialFunction:
					{				
						coupling->setRadialKernel(water_site);
				
						break;
					}
					case ConvCouplingBinaryKernel:
					{				
						coupling->setBinaryKernel(water_site);
						break;
					}
					case ConvCouplingNone: break; //right now unused, but could be useful someday.
					default:
					{
						die("Unknown convolution coupling model specified for water.\n")
					}
				}	
			}
		}
		
		if (params.fluidType != FluidScalarEOSCustom)
		{
		
			//Set the electronic site density model for the coupling
			SiteProperties& Osite = *water.indexedSite[0];
			Osite.siteName="O";
			Osite.couplingZnuc = params.oxygenZnuc;
			Osite.atomicNumber = 8;
						
			SiteProperties& Hsite = *water.indexedSite[1];
			Hsite.siteName="H";
			Hsite.couplingZnuc = params.hydrogenZnuc;
			Hsite.atomicNumber = 1;
		
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
		
		//Create van der Waals mixing functional
		if(params.VDWCouplingScale)
		{	assert(e.vanDerWaals);
			vdwCoupling = new VDWCoupling(*fluidMixture, e.vanDerWaals);
			vdwCoupling->exCorr = &(params.exCorr);
			vdwCoupling->scaleFac = &(params.VDWCouplingScale);
		}
		
		fluidMixture->setPressure(1.01325*Bar);
		
		//---- G=0 constraint -----
		fluidMixture->d0calc = zeroCenter;
	}

	~ConvolutionJDFT()
	{	
		delete coupling;
		delete vdwCoupling;
		delete idgas;
		delete trans;
		delete quad;
		delete fex;
		delete fluidMixture;
		
		if(AwaterCached) delete AwaterCached;
		if(Adiel_rhoExplicitTildeCached) delete Adiel_rhoExplicitTildeCached;
		if(NCached) delete NCached;
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
				ostringstream oss; oss << "N_" << moleculeName;
				sprintf(filename, filenamePattern, oss.str().c_str());
				logPrintf("Dumping %s... ", filename); logFlush();
				saveRawBinary(N[c.offsetDensity], filename); logPrintf("Done.\n"); logFlush();
			}
			else
			{
				for(int j=0; j<c.molecule->nIndices; j++)
				{
					const SiteProperties& s = *c.indexedSite[j];
					ostringstream oss; oss << "N_" << moleculeName << "_" << s.siteName;
					sprintf(filename, filenamePattern, oss.str().c_str());
					logPrintf("Dumping %s... ", filename); logFlush();
					saveRawBinary(N[c.offsetDensity+j], filename); logPrintf("Done.\n"); logFlush();
				}
			}
		}	
	}

	void dumpDebug(const char* filenamePattern) const
	{	
		DataRptrCollection N; char filename[256];
		FluidMixture::Outputs outputs(&N,0,0);		
		fluidMixture->getFreeEnergy(outputs);
		
		//coupling.dumpDebug(filenamePattern);
				
		//loop over sites
		//Compute sphericalized site densities
		for(unsigned ic=0; ic<fluidMixture->get_nComponents(); ic++)
		{
			const FluidMixture::Component& c = fluidMixture->get_component(ic);
			string moleculeName = c.molecule->name;
			if (c.molecule->nIndices == 1)
			{
				ostringstream oss; oss << "Nspherical_" << moleculeName;
				sprintf(filename, filenamePattern, oss.str().c_str());
				logPrintf("Dumping %s... ", filename); logFlush();	
				saveSphericalized(&N[c.offsetDensity], 1, filename);
				logPrintf("Done.\n"); logFlush();
			}
			else
			{
				for(int j=0; j<c.molecule->nIndices; j++)
				{
					const SiteProperties& s = *c.indexedSite[j];
					ostringstream oss; oss << "Nspherical_" << moleculeName << "_" << s.siteName;
					sprintf(filename, filenamePattern, oss.str().c_str());
					logPrintf("Dumping %s... ", filename); logFlush();
					saveSphericalized(&N[c.offsetDensity+j], 1, filename);
					logPrintf("Done.\n"); logFlush();
				}
			}
		}
				
		
		/*
		//Compute potential on O and H sites (done here instead of inside coupling, since need access to densities)
		DataRptrOH N, Adiel_N; fluidMixture->getFreeEnergy(&N); EnergyComponents omega;
		DataGptrOH Ntilde=J(N), Adiel_Ntilde; DataRptrVec p, Adiel_p, eps, Adiel_eps;
		MomentFields mf = {omega, N, Ntilde, p, eps, Adiel_N, Adiel_Ntilde, Adiel_p, Adiel_eps};
		nullToZero(Adiel_N, e.gInfo); nullToZero(Adiel_Ntilde, e.gInfo);
		((ConvolutionJDFT*)this)->coupling.calculateCoupling(mf); //and now Adiel_N and Adiel_Ntilde contain the coupling potential!
		Adiel_N += Jdag(Adiel_Ntilde); Adiel_Ntilde={0,0}; //now all the coupling potential is in real space

		char filename[256];
		sprintf(filename, filenamePattern, "VO");
		logPrintf("Dumping %s... ", filename); logFlush();
		saveRawBinary(Adiel_N.O(), filename); logPrintf("Done.\n"); logFlush();
		sprintf(filename, filenamePattern, "VH");
		logPrintf("Dumping %s... ", filename); logFlush();
		saveRawBinary(Adiel_N.H(), filename); logPrintf("Done.\n"); logFlush();
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
			Adiel_rhoExplicitTildeCached = new DataGptr;
			nullToZero(*Adiel_rhoExplicitTildeCached,e.gInfo);			
			if(vdwCoupling)
			{
				NCached = new DataRptrCollection;				
				nullToZero(*NCached,e.gInfo,fluidMixture->get_nDensities());
			}
			
			FluidMixture::Outputs outputs(NCached,0,Adiel_rhoExplicitTildeCached);
			
			*AwaterCached = fluidMixture->getFreeEnergy(outputs); //Fluid free energy including coupling
			*AwaterCached -= dot(rhoExplicitTilde,O(*Adiel_rhoExplicitTildeCached)); //subtract Electrostatic coupling energy
			*AwaterCached -= coupling->computeElectronic(); //subtract Nonlinear coupling energy
			if(vdwCoupling)
				*AwaterCached -= vdwCoupling->computeElectronic(NCached);
			
		}
		
	}

	void minimizeFluid()
	{	
		assert(AwaterCached); //Ensure that set() was called before calling minimize_fluid()
		TIME("Fluid minimize", globalLog,
			fluidMixture->minimize(e.fluidMinParams);
			
			FluidMixture::Outputs outputs(NCached,0,Adiel_rhoExplicitTildeCached);
			
			*AwaterCached = fluidMixture->getFreeEnergy(outputs); //Fluid free energy including coupling
			*AwaterCached -= dot(fluidMixture->rhoExternal,O(*Adiel_rhoExplicitTildeCached)); //subtract Electrostatic coupling energy
			*AwaterCached -= coupling->computeElectronic(); //subtract Nonlinear coupling energy
			if(vdwCoupling)
				*AwaterCached -= vdwCoupling->computeElectronic(NCached);
		)
	}
	
	double get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
	{
		assert(AwaterCached); //Ensure that set() was called before calling get_Adiel_and_grad()
		
		Adiel_rhoExplicitTilde = clone(*Adiel_rhoExplicitTildeCached);	
		//Uncomment to print the components of Adiel (fluid alone, electrostatic coupling, convolution coupling) 
		
/*		logPrintf("Awater: %le ExtCoulomb: %le ConvCoupling: %le ", *AwaterCached, 
				dot(fluidMixture->rhoExternal,O(*Adiel_rhoExplicitTildeCached)), coupling->computeElectronic());
		if(vdwCoupling)
			logPrintf("VDWCoupling: %le\n", vdwCoupling->computeElectronic(NCached));
		else
			logPrintf("\n"); 
*/		
		double Adiel = *AwaterCached + dot(fluidMixture->rhoExternal,O(*Adiel_rhoExplicitTildeCached))
				+ coupling->computeElectronic(&Adiel_nCavityTilde);
		if(vdwCoupling)
			Adiel += vdwCoupling->computeElectronic(NCached, &extraForces);
		return Adiel;
	}
};



//---------------------------------------------------------------------
//----------------  Interface to the electronic code ------------------
//---------------------------------------------------------------------

FluidSolver::FluidSolver(const Everything& e) : e(e)
{
}

FluidSolver* createFluidSolver(const Everything& e, FluidSolverParams& params)
{	if(params.fluidType != FluidNone)
		Citations::add("Framework of Joint Density Functional Theory", "S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)");
	logPrintf("%s", params.initWarnings.c_str());
	switch(params.fluidType)
	{	case FluidNone:
			return 0; //No solver needed
		case FluidLinearPCM:
			return new LinearPCM(e, params);
		case FluidNonlinearPCM:
			return new NonlinearPCM(e, params);
		case FluidNonlocalPCM:
			return new NonlocalPCM(e, params);
		default: //All explicit fluid functionals:
			return new ConvolutionJDFT(e, params);
	}
}
