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

#include <fluid/Fex_H2O_ScalarEOS_internal.h>
#include <fluid/Fex_H2O_ScalarEOS.h>
#include <fluid/Fex_LJ.h>
#include <core/Units.h>
#include <core/Operators.h>

string rigidMoleculeCDFT_ScalarEOSpaper =  "R. Sundararaman and T.A. Arias, (to be submitted to Comp. Phys. Comm.)";

static const double rOH = 1.0*Angstrom; 
static const double thetaHOH = acos(-1.0/3); 
//static const double rOH = 0.96719*Angstrom; //HACK WARNING
//static const double thetaHOH = 1.8069; //103.53 degrees HACK WARNING

Fex_H2O_ScalarEOS::Fex_H2O_ScalarEOS(FluidMixture& fluidMixture)
: Fex(fluidMixture),
fex_LJatt(gInfo), siteChargeKernel(gInfo),
eval(new ScalarEOS_eval(T)),
propO(gInfo, eval->sphereRadius,0.0, 0.8476,&siteChargeKernel),
propH(gInfo, 0.0,0.0,-0.4238,&siteChargeKernel),
molecule("H2O",
	&propO,
		 vector3<>(0,0,0),
	&propH,
		 vector3<>(0, -rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)),
		 vector3<>(0, +rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)) )
{
	//Initialize the kernels:
	applyFuncGsq(gInfo, setCoulombCutoffKernel, siteChargeKernel.data); siteChargeKernel.set();
	setLJatt(fex_LJatt, -9.0/(32*sqrt(2)*M_PI*pow(2*eval->sphereRadius,3)), 2*eval->sphereRadius);
	Citations::add("Scalar-EOS water functional", rigidMoleculeCDFT_ScalarEOSpaper);
}

double Fex_H2O_ScalarEOS::get_aDiel() const
{
	return 1 - T/(7.35e3*Kelvin); 
}


#ifdef GPU_ENABLED
void Fex_H20_ScalarEOS_gpu(int nr, const double* Nbar, double* Fex, double* grad_Nbar, ScalarEOS_eval eval);
#endif
double Fex_H2O_ScalarEOS::compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const
{	//Compute LJatt weighted density:
	DataRptr Nbar = I(Ntilde[0]*fex_LJatt), grad_Nbar; nullToZero(grad_Nbar, gInfo);
	//Evaluated weighted density functional:
	DataRptr Aex(DataR::alloc(gInfo,isGpuEnabled()));
	#ifdef GPU_ENABLED
	Fex_H20_ScalarEOS_gpu(gInfo.nr, Nbar->dataGpu(), Aex->dataGpu(), grad_Nbar->dataGpu(), *eval);
	#else
	threadedLoop(eval, gInfo.nr, Nbar->data(), Aex->data(), grad_Nbar->data());
	#endif
	//Convert gradients:
	DataRptr NO = I(Ntilde[0]);
	grad_Ntilde[0] += fex_LJatt*Idag(NO*grad_Nbar) + Idag(Aex);
	return gInfo.dV*dot(NO,Aex);
}

double Fex_H2O_ScalarEOS::computeUniform(const double* N, double* grad_N) const
{	double AexPrime, Aex;
	(*eval)(0, &N[0], &Aex, &AexPrime);
	grad_N[0] += Aex + N[0]*AexPrime;	
	return N[0]*Aex;
}

std::vector<std::vector<vector3<>>> getPositionList(std::vector<H2OSite>& H2OSites)
{
	int OsiteCounter = 0, HsiteCounter = 0, siteCounter = 0, siteIndex = 0;	
	std::vector<std::vector<vector3<>>> PositionList;
	PositionList.resize(H2OSites.size(),H2OSites[0].Positions);
	
	for (uint iSite=0; iSite<H2OSites.size(); iSite++)
	{	
		if(H2OSites[iSite].name == "O")
		{
			//assumes oxygen is first in the list of site properties
			siteIndex = 0;
			OsiteCounter = H2OSites[iSite].Positions.size();
		}
		else if(H2OSites[iSite].name == "H")
		{
			//assumes hydrogen is second in the list of site properties
			siteIndex = 1;
			HsiteCounter =  H2OSites[iSite].Positions.size();
		}
		else
		{
			siteCounter++;
			siteIndex = siteCounter+1;
		}
		
		PositionList[siteIndex] = H2OSites[iSite].Positions;
		
		logPrintf("Initialized water sites of type %s at positions (in Bohr):\n", H2OSites[iSite].name.c_str());
		for (uint iPos=0; iPos<PositionList[siteIndex].size(); iPos++)
			logPrintf("\t %lg %lg %lg \n",PositionList[siteIndex][iPos][0] ,PositionList[siteIndex][iPos][1], PositionList[siteIndex][iPos][2]); 
	}
	
	if (!((HsiteCounter == 2) && (OsiteCounter==1)))
		die("Need to specify one O site and two H sites for custom designed ScalarEOS water functional.");
	
	return PositionList;
} 

std::vector<SiteProperties*> getSitePropList(const GridInfo& gInfo, std::vector<H2OSite>& H2OSites, RealKernel* siteChargeKernelPtr, const double& OsphereRadius)
{
	bool Osite = false, Hsite = false;
	int siteCounter = 0,siteIndex = 0;
	double sphereRadius = 0.0;
	double sphereSigma=0.0;
	std::vector<SiteProperties*> PropList;
	PropList.resize(H2OSites.size(),0);
	
	for (uint iSite=0; iSite<H2OSites.size(); iSite++)
	{
		if(H2OSites[iSite].name == "O")
		{
			//assumes oxygen is first in the list of site properties
			siteIndex = 0;
			sphereRadius = OsphereRadius;
			Osite=true;
		}
		else if(H2OSites[iSite].name == "H")
		{
			//assumes hydrogen is second in the list of site properties
			siteIndex = 1;
			sphereRadius = 0.0;
			Hsite=true;
		}
		else
		{
			siteCounter++;
			sphereRadius = 0.0;
			siteIndex = siteCounter+1;
		}
		
		logPrintf("Custom ScalarEOS water site %s initialized with sphere radius %lg and sigma %lg\n", H2OSites[iSite].name.c_str(), sphereRadius, sphereSigma);
		
					
		PropList[siteIndex] = new SiteProperties(gInfo, sphereRadius, sphereSigma, H2OSites[iSite], siteChargeKernelPtr);
	}	
	
	if (!(Hsite && Osite))
		die("Need to specify both O and H sites for custom designed ScalarEOS water functional");
		
	return PropList;
} 

double calc_aDielFactor(double dipoleMoment)
{
	//use function fit from octave for dipole moments between 0 and 1.5 ea_0
	double polyCoeff[] = {4.4266e+04,  -2.9262e+05,   8.3304e+05,  -1.3257e+06,   1.2890e+06,  -7.8496e+05,
	2.9834e+05,  -6.6193e+04,   1.4463e+04,  -4.1945e+02,   6.1205e+00};
	int nCoeff = 11;
	double factor=0.0;
	for (int i=0; i<nCoeff; i++)
	{
		factor += polyCoeff[nCoeff-i-1]*pow(dipoleMoment,double(i));
	}
	return factor;
}

Fex_H2O_Custom::Fex_H2O_Custom(FluidMixture& fluidMixture, std::vector<H2OSite>& H2OSites)
: Fex_H2O_ScalarEOS(fluidMixture), prop(getSitePropList(gInfo, H2OSites, &siteChargeKernel, eval->sphereRadius)),
pos(getPositionList(H2OSites)),customMolecule(prop, pos, "H2O")
{
			//put remainder of constructor here
			if ( fabs(customMolecule.get_charge()) > 1e-14)
				die("Custom designed ScalarEOS water functional requires neutral water molecule");
			
			double dipoleMoment = fabs(customMolecule.get_dipole());
			
			aDielFactor = calc_aDielFactor(dipoleMoment);
			
			logPrintf("Custom ScalarEOS functional created with charge: %lg, dipole moment: %lg and aDiel factor: %lg\n",
					  customMolecule.get_charge(),dipoleMoment,aDielFactor);
				
}

double Fex_H2O_Custom::get_aDiel() const
{
	return 1 - T/(aDielFactor*Kelvin); 
}

