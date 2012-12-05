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

#include <electronic/ConvCoupling.h>
#include <electronic/ExCorr.h>
#include <fluid/FluidMixture.h>
#include <fluid/Molecule.h>
#include <core/DataIO.h>
#include <core/Data.h>
#include <electronic/RadialFunction.h>
#include <electronic/operators.h>

ConvCoupling::ConvCoupling(FluidMixture& fluidMixture)
: Fmix(fluidMixture),CouplingData(gInfo)
{
	//initialize nFluid
	//nullToZero(nFluid,gInfo);
	//initialize the pointer in FluidMixture to point to the CouplingData structure in ConvCoupling
	fluidMixture.ConvCouplingPtr = &CouplingData;
	Citations::add("Convolution-coupling for Joint Density Functional Theory",
		"K. Letchworth-Weaver, R. Sundararaman and T.A. Arias, (under preparation)");
}


ConvCoupling::~ConvCoupling()
{
	//loop over all sites and deallocate the coupling kernels
	for(unsigned ic=0; ic<fluidMixture.get_nComponents(); ic++)
	{		
		const FluidMixture::Component& c = fluidMixture.get_component(ic);
		for(int j=0; j<c.molecule->nIndices; j++)
		{
			const SiteProperties& s = *c.indexedSite[j];
			delete s.couplingElecKernel;			
		}
	}
	
}

double ConvCoupling::computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const
{	return 0.0; //No electronic systen to couple to in the bulk fluid
}

double ConvCoupling::compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const 
{	
	
	/*DataGptr nFluidTilde;
	nullToZero(nFluidTilde,gInfo);
		for(unsigned ic=0; ic<fluidMixture.get_nComponents(); ic++)
	{
		const FluidMixture::Component& c = fluidMixture.get_component(ic);
		for(int j=0; j<c.molecule->nIndices; j++)
		{
			const SiteProperties& s = *c.indexedSite[j];
			if(s.couplingZnuc && s.couplingElecKernel)
			{
				nFluidTilde += (*s.couplingElecKernel * Ntilde[c.offsetDensity+j]);
			}		
		}
	}*/

	
	if(CouplingData.nFluidTilde==0)
	{
		DataGptr& tmp = ((ConvCoupling*)this)->CouplingData.nFluidTilde; //caching data, so casting non-const
		nullToZero(tmp,gInfo);
		logPrintf("WARNING: nFluidTilde not allocated in ConvCoupling.compute()!\n"); 
	}
	
	
	/*DataRptr& nFluid = ((ConvCoupling*)this)->nFluid; //caching data, so casting non-const
	nFluid = I(nFluidTilde);*/
	
	
	//NOTE: Assumes that nFluidTilde has already been calculated by operator() in FluidMixture
	DataRptr nFluid = I(CouplingData.nFluidTilde);
	
			
	//declare energy and gradient variables
	double PhiCoupling;
	DataGptr grad_nFluidTilde;
	DataRptr nTot = nFluid + nCavity;

	//Calculate exchange, correlation, and kinetic energy
	DataRptr Vxc_tot, Vxc_fluid;
	PhiCoupling =
	+ (*exCorr)(nTot, &Vxc_tot, true)
	- (*exCorr)(nFluid, &Vxc_fluid, true)
	- (*exCorr)(nCavity, 0,  true);
	//PhiCoupling *= CouplingScale;
	
	grad_nFluidTilde = O(J(Vxc_tot - Vxc_fluid));
	//grad_nFluidTilde = CouplingScale*O(J(Vxc_tot - Vxc_fluid));
	
	
	//loop over sites and use the chain rule to calculate nonlinear contribution to the coupling.
	for(unsigned ic=0; ic<fluidMixture.get_nComponents(); ic++)
	{		
		const FluidMixture::Component& c = fluidMixture.get_component(ic);
		for(int j=0; j<c.molecule->nIndices; j++)
		{
			const SiteProperties& s = *c.indexedSite[j];
			if(s.couplingElecKernel)
			{
				grad_Ntilde[c.offsetDensity+j] += 1.0/gInfo.dV*(*s.couplingElecKernel *grad_nFluidTilde);
			}	
		}
	}
	
	return PhiCoupling;
	}

string ConvCoupling::getName() const
{	return "ConvCoupling";
}

//function for setting an exponential kernel for a given site with charge Z and width a => Z/(8*pi*a^3)*exp(-r/a)
inline void setExpKernelFunc(int i, double G2, double* kernel, double Z, double a)
{	
	kernel[i] = Z/pow(1 + a*a*G2, 2);
}

//function for setting an exponential kernel for a given site with charge Z and width a => Z/(8*pi*a^3)*exp(-r/a)
inline void setExpCusplessKernelFunc(int i, double G2, double* kernel, double Z, double a)
{	
	kernel[i] = Z/pow(1 + a*a*G2, 3);
}

void ConvCoupling::setExponentialKernel(SiteProperties& s)
{	
	//allocate the kernel
	s.couplingElecKernel = new RealKernel(gInfo);
	
	double elecCharge = s.couplingZnuc+s.convCouplingSiteCharge;
	logPrintf("Created exponential electron density model for %s site with width %.12lf and norm %.12lf.\n",
			  s.siteName.c_str(),s.convCouplingWidth,elecCharge);
	if (s.kernelFilename.length() == 0) //If there is not another filename to be added in addition to this one.
	{
		if (fabs(s.convCouplingSiteCharge-s.chargeZ*s.chargeKernel->data[0])>1e-12)
			logPrintf("Warning: classical site charges and electron kernel charges not balanced in convolution coupling.\n");
	}
	applyFuncGsq(gInfo, setExpKernelFunc, s.couplingElecKernel->data,elecCharge,s.convCouplingWidth);
	s.couplingElecKernel->set();
}

void ConvCoupling::setExpCusplessKernel(SiteProperties& s)
{	
	//allocate the kernel
	s.couplingElecKernel = new RealKernel(gInfo);
	
	double elecCharge = s.couplingZnuc+s.convCouplingSiteCharge;
	logPrintf("Created cuspless exponential electron density model for %s site with width %.12lf and norm %.12lf.\n",
			  s.siteName.c_str(),s.convCouplingWidth,elecCharge);
	if (s.kernelFilename.length() == 0) //If there is not another filename to be added in addition to this one.
	{
		if (fabs(s.convCouplingSiteCharge-s.chargeZ*s.chargeKernel->data[0])>1e-12)
			logPrintf("Warning: classical site charges and electron kernel charges not balanced in convolution coupling.\n");
	}
	applyFuncGsq(gInfo, setExpCusplessKernelFunc, s.couplingElecKernel->data,elecCharge,s.convCouplingWidth);
	s.couplingElecKernel->set();
}

void ConvCoupling::setBinaryKernel(SiteProperties& s)
{
	//allocate the kernel
	s.couplingElecKernel = new RealKernel(gInfo);
	
	//Read kernel in full G space
	complexDataGptr Kernel_fullG(complexDataG::alloc(gInfo));
	loadRawBinary(Kernel_fullG, s.kernelFilename.c_str());
	//reduce to half G space:
	DataGptr Kernel = Real(Kernel_fullG);
	double Kernel0 = Kernel->data()[0].real();
	logPrintf("Read electron density model for %s site from binary file '%s', with net charge %.12lf and site charge %.12lf\n",
				s.siteName.c_str(), s.kernelFilename.c_str(), Kernel0, Kernel0-s.couplingZnuc);
	const complex* KernelData = Kernel->data();
	if (fabs(Kernel0-s.couplingZnuc-s.chargeZ*s.chargeKernel->data[0])>1e-12)
		logPrintf("Warning: classical site charges and electron kernel charges not balanced in convolution coupling.\n");
	for(int i=0; i<gInfo.nG; i++)
	{
		s.couplingElecKernel->data[i] = KernelData[i].real();
	}
	s.couplingElecKernel->set();	
}

void ConvCoupling::setRadialKernel(SiteProperties& s)
{	
	ifstream ifs(s.kernelFilename.c_str());
	if(!ifs.is_open()) die("Can't open radial electron density file '%s' for reading.\n", s.kernelFilename.c_str());
	logPrintf("\nReading radial electron density model for %s site from '%s':\n",s.siteName.c_str(),s.kernelFilename.c_str());
	 
	std::vector<double> rVec,nVec; //to store values of radius and density
	double deltaRMin = 10.0; //minimum distance between radial grid points.
	
	while (!ifs.eof())
	{
		double r,n;
		string line;
		getline(ifs, line);
		istringstream(line) >> r >> n;
		logPrintf("r: %.12lf n: %.12lf",r,n);
		rVec.push_back(r);
		nVec.push_back(n);	
		
		int Nr = rVec.size();
		double deltaR = rVec[Nr]-rVec[Nr-1];
		
		if (deltaR <= 0.0)
		{
			die("ERROR reading electron density model for %s site from %s:\n"
				"Radial gridpoints must be in ascending order\n", s.siteName.c_str(), s.kernelFilename.c_str());
		}
		else if (deltaR < deltaRMin)
			deltaRMin = deltaR;
		
		if( ifs.fail() || ifs.bad())
			die("ERROR reading electron density model for %s site from %s\n", s.siteName.c_str(), s.kernelFilename.c_str());
		
	}
	
	RadialFunctionR KernelR(rVec.size());
	KernelR.r = rVec;
	KernelR.f = nVec;
	KernelR.initWeights();
	
	double dG = 0.02; //The uniform G-spacing for the radial function
	int ngridG = 2 * ceil(2.0*M_PI/(deltaRMin*dG)); //set the number of G grid points
	RadialFunctionG KernelG;
	
	//transform to fourier space
	KernelR.transform(0, dG, ngridG, KernelG);
	
	//allocate the kernel if not already allocated
	if(!s.couplingElecKernel)
	{
		s.couplingElecKernel = new RealKernel(gInfo);	
		//Interpolate fourier space function onto gInfo
		radialFunctionG(KernelG, *s.couplingElecKernel);
	}
	else //assumes radial function is an addition to an exponential (or other analytical function) 
	{
		RealKernel tmp(gInfo);	
		//Interpolate fourier space function onto gInfo
		radialFunctionG(KernelG, tmp);
		for(int i=0; i<gInfo.nG; i++)
			s.couplingElecKernel->data[i] += tmp.data[i];
	}
	
	double Kernel0 = s.couplingElecKernel->data[0];
	if (fabs(Kernel0-s.couplingZnuc-s.chargeZ*s.chargeKernel->data[0])>1e-12)
		logPrintf("Warning: classical site charges and electron kernel charges not balanced in convolution coupling.\n");

	logPrintf("Final electron density model for %s site has net charge %.12lf and site charge %.12lf\n",
				s.siteName.c_str(), Kernel0, Kernel0-s.couplingZnuc);
				
	s.couplingElecKernel->set();	
}

void ConvCoupling::setExplicit(const DataRptr& nCavity)
{
	this->nCavity = nCavity;
}

void ConvCoupling::setExplicit(const DataGptr& nCavityTilde)
{
	this->nCavity = clone(I(nCavityTilde));	
}

double ConvCoupling::computeElectronic(DataGptr* grad_nCavityTilde)
{
	//only returns the nonlinear pieces
	
	if(CouplingData.nFluidTilde==0)
	{
		DataGptr& tmp = ((ConvCoupling*)this)->CouplingData.nFluidTilde; //caching data, so casting non-const
		nullToZero(tmp,gInfo);
		logPrintf("WARNING: nFluidTilde not allocated in ConvCoupling.computeElectronic()!\n");
	}
			
	DataRptr nFluid = I(CouplingData.nFluidTilde);
	
	/*if(!nFluid)
	{
		DataRptr& tmp = ((ConvCoupling*)this)->nFluid; //caching data, so casting non-const
		nullToZero(tmp,gInfo);
		logPrintf("WARNING: nFluid not allocated in ConvCoupling.computeElectronic()!\n");
	}*/
	
	DataRptr nTot = nFluid + nCavity;
	DataRptr Vxc_tot, Vxc_cavity;
	
	double Acoupling = 
		+ (*exCorr)(nTot, grad_nCavityTilde ? &Vxc_tot : 0, true)
		- (*exCorr)(nCavity, grad_nCavityTilde ? &Vxc_cavity : 0, true)
		- (*exCorr)(nFluid, 0, true);
		
		if(grad_nCavityTilde)
		*grad_nCavityTilde = J(Vxc_tot - Vxc_cavity);

	return Acoupling;
}

void ConvCoupling::dumpDebug(const char* filenamePattern) const
{
	
}

