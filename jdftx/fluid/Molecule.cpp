/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <fluid/ErfFMTweight.h>
#include <fluid/Molecule.h>
#include <core/Operators.h>
#include <electronic/ColumnBundle.h>


Molecule::Site::Site(string name, int atomicNumber) : name(name), Rhs(0), atomicNumber(atomicNumber), Znuc(0), sigmaNuc(0), Zelec(0), aElec(0), alpha(0), aPol(0), initialized(false)
{	
}

Molecule::Site::~Site()
{	free();
}

void Molecule::Site::free()
{	if(initialized)
	{	elecKernel.free();
		chargeKernel.free();
		polKernel.free();
		if(Rhs)
		{	w0.free();
			w1.free();
			w2.free();
			w3.free();
			w1v.free();
			w2m.free();
		}
	}
}

//Fourier transform of cuspless exponential
inline double cusplessExpTilde(double G, double norm, double a)
{	double aG = a*G;
	double den = 1./(1.+aG*aG);
	return norm * den*den*den;
}

//Fourier transform of gaussian
inline double gaussTilde(double G, double norm, double sigma)
{	double sigmaG = sigma*G;
	return norm * exp(-0.5*sigmaG*sigmaG);
}

void Molecule::Site::setup(const GridInfo& gInfo)
{	if(initialized) free();
	logPrintf("     Initializing site '%s'\n", name.c_str());
	
	//Initialize electron density kernel:
	if(elecFilename.length() || Zelec)
	{	logPrintf("       Electron density: ");
		if(elecFilename.length())
		{	logPrintf("reading from '%s' ...", elecFilename.c_str());
			//TODO: read from file
		}
		else
		{	logPrintf("cuspless exponential with width %lg and", aElec);
			elecKernel.init(0, gInfo.dGradial, gInfo.GmaxGrid, cusplessExpTilde, Zelec, aElec);
		}
		logPrintf(" norm %lg\n", elecKernel(0));
	}
	
	//Initialize charge kernel:
	if(elecKernel || Znuc)
	{	logPrintf("       Charge density: gaussian nuclear width %lg", sigmaNuc);
		std::vector<double> samples(unsigned(ceil(gInfo.GmaxGrid/gInfo.dGradial))+5, 0.);
		for(unsigned iG=0; iG<samples.size(); iG++)
		{	double G = iG * gInfo.dGradial;
			if(elecKernel) samples[iG] += elecKernel(G);
			if(Znuc) samples[iG] -= gaussTilde(G, Znuc, sigmaNuc);
		}
		chargeKernel.init(0, samples, gInfo.dGradial);
		logPrintf(" with net site charge %lg\n", chargeKernel(0));
	}
	
	//Initialize polarizability kernel:
	if(alpha)
	{	logPrintf("       Polarizability: cuspless exponential with width %lg and norm %lg\n", aPol, alpha);
		polKernel.init(0, gInfo.dGradial, gInfo.GmaxGrid, cusplessExpTilde, 1., aPol);
	}
	
	if(Rhs)
	{	logPrintf("       Hard sphere radius: %lg bohrs\n", Rhs);
		//Initialize hard sphere weights:
		ErfFMTweight erfFMTweight(Rhs, 0.);
		unsigned nGradial = unsigned(ceil(gInfo.GmaxGrid/gInfo.dGradial)+5);
		std::vector<double> w0(nGradial), w1(nGradial), w2(nGradial), w3(nGradial), w1v(nGradial), w2m(nGradial);
		for(unsigned iG=0; iG<nGradial; iG++)
			erfFMTweight(iG*gInfo.dGradial, w0[iG], w1[iG], w2[iG], w3[iG], w1v[iG], w2m[iG]);
		this->w0.init(0, w0, gInfo.dGradial);
		this->w1.init(0, w1, gInfo.dGradial);
		this->w2.init(0, w2, gInfo.dGradial);
		this->w3.init(0, w3, gInfo.dGradial);
		this->w1v.init(0, w1v, gInfo.dGradial);
		this->w2m.init(0, w2m, gInfo.dGradial);
	}
	
	logPrintf("       Positions in reference frame:\n");
	for(vector3<> r: positions) logPrintf("         [ %+.6lf %+.6lf %+.6lf ]\n", r[0], r[1], r[2]);
	initialized = true;
}


/*
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
*/

Molecule::Molecule(string name) : name(name), initialized(false)
{
}

Molecule::~Molecule()
{	if(initialized)
	{	mfKernel.free();
	}
}

inline double sphericalShellTilde(double G, double R)
{	return bessel_jl(0, G*R);
}

void Molecule::setup(const GridInfo& gInfo, double Rmf)
{	logPrintf("   Initializing fluid molecule '%s'\n", name.c_str());
	for(auto& site: sites) site->setup(gInfo);
	logPrintf("     Net charge: %lg   dipole magnitude: %lg\n", getCharge(), getDipole().length());
	mfKernel.init(0, gInfo.dGradial, gInfo.GmaxGrid, sphericalShellTilde, Rmf ? Rmf : pow(3*getVhs()/(4.*M_PI), 1./3));
	initialized = true;
}


bool Molecule::isMonoatomic() const
{	return (sites.size()==1) && (sites[0]->positions.size()==1);
}

double Molecule::getCharge() const
{	double Q = 0.0;
	for(const auto& site: sites)
		if(site->chargeKernel)
			Q += site->chargeKernel(0) * site->positions.size();
	if(fabs(Q) < 1e-12) return 0.; //simplify neutrality testing
	else return Q;
}

vector3<> Molecule::getDipole() const
{	vector3<> P;
	for(const auto& site: sites)
		if(site->chargeKernel)
			for(const vector3<>& r: site->positions)
				P += site->chargeKernel(0) * r;
	if(P.length() < 1e-12) return vector3<>(); //simplify polarity testing
	else return P;
}

double Molecule::getVhs() const
{	double Vhs = 0;
	for(const auto& site: sites)
		if(site->w3)
			Vhs += site->w3(0) * site->positions.size();
	return Vhs;
}

double Molecule::getAlphaTot() const
{	double alphaTot = 0.;
	for(const auto& site: sites)
		if(site->polKernel)
			alphaTot += site->alpha * site->positions.size();
	return alphaTot;
}

std::map<double,int> Molecule::getBonds() const
{	std::map<double,int> bond;
	for(const auto& site1: sites)
	{	double R1 = site1->Rhs;
		if(R1) for(vector3<> pos1: site1->positions)
		{	for(const auto& site2: sites)
			{	double R2 = site2->Rhs;
				if(R2) for(vector3<> pos2: site2->positions)
				{	if(fabs(R1+R2-(pos1-pos2).length()) < 1e-6*(R1+R2))
						bond[R1*R2/(R1+R2)]++;
				}
			}
		}
	}
	//correct for double counting:
	for(auto& bondEntry: bond)
	{	assert(bondEntry.second % 2 == 0);
		bondEntry.second /= 2;
	}
	return bond;
}

void Molecule::setModelMonoatomic(string name, double Q, double Rhs)
{	sites.clear();
	this->name = name;
	auto site = std::make_shared<Molecule::Site>(name);
		site->Znuc = Q; site->sigmaNuc = (1./6)*Rhs;
		site->Rhs = Rhs;
	sites.push_back(site);
	site->positions.push_back(vector3<>());
}
