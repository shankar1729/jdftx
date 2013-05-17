/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Kendra Letchworth Weaver

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

#include <electronic/VanDerWaals.h>
#include <electronic/Everything.h>
#include <electronic/SpeciesInfo_internal.h>
#include <electronic/operators.h>
#include <core/DataMultiplet.h>
#include <core/Units.h>

const static int atomicNumberMaxGrimme = 54;
const static int atomicNumberMax = 118;

//vdW correction energy upto a factor of -s6 (where s6 is the ExCorr dependnet scale)
//for a pair of atoms separated by r, given the C6 and R0 parameters for pair.
//The corresponding derivative w.r.t r is stored in E_r
inline double vdwPairEnergyAndGrad(double r, double C6, double R0, double& E_r)
{
	//Regularize the spurious r=0 singularity in the Grimme vdw functional.
	double invR0 = 1./R0, rByR0 = r * invR0;
	if(rByR0 < 0.3000002494598603)
	{	E_r = 0.;
		return 0.00114064201325433 * C6 * pow(invR0,6);
	}
	
	const double d = 20.;
	double exponential = exp(-d*(rByR0-1.));
	double fdamp = 1./(1. + exponential);
	double fdamp_r = (fdamp*fdamp) * exponential * (d*invR0);
	
	double invr = 1./r;
	double C6invr6 = C6 * pow(invr, 6);
	double C6invr6_r = (-6./r) * C6invr6;
	
	E_r = C6invr6_r * fdamp + C6invr6 * fdamp_r;
	return C6invr6 * fdamp;
}

VanDerWaals::VanDerWaals(const Everything& everything)
{
	logPrintf("\nInitializing van der Waals corrections ... ");
	e = &everything;
	
	// Checks whether pseudopotentials contain atomic numbers
	for(auto sp: e->iInfo.species)
	{	assert(sp->atomicNumber);
		if(sp->atomicNumber > atomicNumberMax) die("\nAtomic numbers > %i not supported!\n", atomicNumberMax);
		if(sp->atomicNumber > atomicNumberMaxGrimme)
			logPrintf("WARNING: Using extended vdW params (beyond Grimme's data set) for species %s with atomic number > %d.\n",
				sp->name.c_str(), atomicNumberMaxGrimme);
	}

	// Constructs the EXCorr -> scaling factor map
	scalingFactor["gga-PBE"] = 0.75;
	scalingFactor["hyb-gga-xc-b3lyp"] = 1.05;
	scalingFactor["mgga-TPSS"] = 1.;

	// Sets up the C6 and R0 parameters
	atomParams.resize(atomicNumberMax+1);
	atomParams[1] = AtomParams(0.14 , 1.001 );
	atomParams[2] = AtomParams(0.08 , 1.012 );
	atomParams[3] = AtomParams(1.61 , 0.825 );
	atomParams[4] = AtomParams(1.61 , 1.408 );
	atomParams[5] = AtomParams(3.13 , 1.485 );
	atomParams[6] = AtomParams(1.75 , 1.452 );
	atomParams[7] = AtomParams(1.23 , 1.397 );
	atomParams[8] = AtomParams(0.70 , 1.342 );
	atomParams[9] = AtomParams(0.75 , 1.287 );
	atomParams[10] = AtomParams(0.63 , 1.243 );
	atomParams[11] = AtomParams(5.71 , 1.144 );
	atomParams[12] = AtomParams(5.71 , 1.364 );
	atomParams[13] = AtomParams(10.79,1.639  );
	atomParams[14] = AtomParams(9.23 , 1.716 );
	atomParams[15] = AtomParams(7.84 , 1.705 );
	atomParams[16] = AtomParams(5.57 , 1.683 );
	atomParams[17] = AtomParams(5.07 , 1.639 );
	atomParams[18] = AtomParams(4.61 , 1.595 );
	atomParams[19] = AtomParams(10.80 , 1.485);
	atomParams[20] = AtomParams(10.80 , 1.474);
	for(int Z=21; Z<=30; Z++)
		atomParams[Z] = AtomParams(10.80 , 1.562);
	atomParams[31] = AtomParams(16.99 , 1.650);
	atomParams[32] = AtomParams(17.10 , 1.727);
	atomParams[33] = AtomParams(16.37 , 1.760);
	atomParams[34] = AtomParams(12.64 , 1.771);
	atomParams[35] = AtomParams(12.47 , 1.749);
	atomParams[36] = AtomParams(12.01 , 1.727);
	atomParams[37] = AtomParams(24.67 , 1.628);
	atomParams[38] = AtomParams(24.67 , 1.606);
	for(int Z=39; Z<=48; Z++)
		atomParams[Z] = AtomParams(24.67 , 1.639);
	atomParams[49] = AtomParams(37.32 , 1.672);
	atomParams[50] = AtomParams(38.71 , 1.804);
	atomParams[51] = AtomParams(38.44 , 1.881);
	atomParams[52] = AtomParams(31.74 , 1.892);
	atomParams[53] = AtomParams(31.50 , 1.892);
	atomParams[54] = AtomParams(29.99 , 1.881);
	//------ Extending Grimme's set to all remaining elements using his formula:
	//(with IP data from Hohm et al., J. Phys. Chem. A 116, 697 (2012))
	atomParams[55] = AtomParams(47 , 1.798);
	atomParams[56] = AtomParams(47 , 1.764);
	for(int Z=57; Z<=80; Z++)
		atomParams[Z] = AtomParams(47 , 1.76);
	atomParams[81] = AtomParams(63.7 , 1.754);
	atomParams[82] = AtomParams(50.8 , 1.816);
	atomParams[83] = AtomParams(51.5 , 1.881);
	atomParams[84] = AtomParams(53.4 , 1.898);
	atomParams[85] = AtomParams(56.3 , 1.915);
	atomParams[86] = AtomParams(54.8 , 1.907);
	atomParams[87] = AtomParams(55 , 1.851);
	atomParams[88] = AtomParams(55 , 1.844);
	for(int Z=89; Z<=atomicNumberMax; Z++)
		atomParams[Z] = AtomParams(55 , 1.75);
	logPrintf("Done!\n");
	
	Citations::add("Van der Waals correction pair-potentials", "S. Grimme, J. Comput. Chem. 27, 1787 (2006)");
}

double VanDerWaals::energyAndGrad(std::vector<Atom>& atoms, const double scaleFac) const
{
	//Truncate summation at 1/r^6 < 10^-16 => r ~ 100 bohrs
	vector3<bool> isTruncated = e->coulombParams.isTruncated();
	vector3<int> n;
	for(int k=0; k<3; k++)
		n[k] = isTruncated[k] ? 0 : (int)ceil(100. / e->gInfo.R.column(k).length());
	
	double Etot = 0.;  //Total VDW Energy
	for(int c1 = 0; c1 < (int) atoms.size(); c1++)
	{	const AtomParams& c1params = getParams(atoms[c1].atomicNumber);
		for(int c2 = 0; c2 < (int) atoms.size(); c2++)
		{	const AtomParams& c2params = getParams(atoms[c2].atomicNumber);
			
			if(c1 == c2) continue; // No self interaction		
			double C6 = sqrt(c1params.C6 * c2params.C6);
			double R0 = c1params.R0 + c2params.R0;
			
			vector3<int> t;
			for(t[0] = -n[0]; t[0]<=n[0]; t[0]++)
			for(t[1] = -n[1]; t[1]<=n[1]; t[1]++)
			for(t[2] = -n[2]; t[2]<=n[2]; t[2]++)
			{
				vector3<> RijVec = e->gInfo.R * (atoms[c2].pos + t - atoms[c1].pos);
				double Rij = RijVec.length();
				vector3<> RijHat = RijVec / Rij;
				
				double E_Rij, E = vdwPairEnergyAndGrad(Rij, C6, R0, E_Rij);
				
				Etot -= 0.5 * scaleFac * E;
				atoms[c1].force += scaleFac * E_Rij * RijHat;
			}
		}
	}
	return Etot;
}


double VanDerWaals::energyAndGrad(const DataGptrCollection& Ntilde, const std::vector< int >& atomicNumber, const double scaleFac,
	DataGptrCollection* grad_Ntilde, IonicGradient* forces) const
{		
	
	double Etot = 0.;
	const GridInfo& gInfo = e->gInfo;
	const std::vector< std::shared_ptr<SpeciesInfo> >& species = e->iInfo.species;
	
	for(unsigned i=0; i<species.size(); i++) //Loop over species of explicit system
	{	
		std::shared_ptr<SpeciesInfo> sp = species[i];
		DataGptr SG(DataG::alloc(gInfo, isGpuEnabled()));
		DataGptr ccgrad_SG; //set grad wrt structure factor
		int nAtoms = sp->atpos.size(); //number of atoms of ith species
		
		callPref(getSG)(gInfo.S, nAtoms, sp->atposPref, 1.0/gInfo.detR, SG->dataPref()); //get structure factor SG for atom type i

		for(unsigned j=0; j<atomicNumber.size(); j++) //Loop over sites in the fluid
			if(atomicNumber[j]) //Check to make sure fluid site should include van der Waals corrections
			{
				const RadialFunctionG& Kernel_ij = getRadialFunction(sp->atomicNumber,atomicNumber[j]); //get ij radial function
				DataGptr E_Ntilde = (-scaleFac * gInfo.nr) * (Kernel_ij * SG); //calculate effect of ith explicit atom on gradient wrt jth site density
				Etot += gInfo.dV * dot(Ntilde[j], E_Ntilde); //accumulate into total energy
				if(grad_Ntilde)
					(*grad_Ntilde)[j] += E_Ntilde; //accumulate into gradient wrt jth site density
				if(forces)
					ccgrad_SG += (-scaleFac) * (Kernel_ij * Ntilde[j]); //accumulate forces on ith atom type from jth site density
			}

		if(forces && ccgrad_SG) //calculate forces due to ith atom
		{	DataGptrVec gradAtpos; nullToZero(gradAtpos, gInfo);
			vector3<complex*> gradAtposData; for(int k=0; k<3; k++) gradAtposData[k] = gradAtpos[k]->dataPref();
			for(int at=0; at<nAtoms; at++)
			{	
				callPref(gradSGtoAtpos)(gInfo.S, sp->atpos[at], ccgrad_SG->dataPref(), gradAtposData);
				for(int k=0; k<3; k++)
					(*forces)[i][at][k] -= sum(gradAtpos[k]); //negative gradient on ith atom type
			}
		}
	}
	
	return Etot;
}

double VanDerWaals::getScaleFactor(string exCorrName, double scaleOverride) const
{	if(scaleOverride) return scaleOverride;
	auto iter = scalingFactor.find(exCorrName);
	if(iter == scalingFactor.end())
		die("\nGrimme vdW scale factor not known for functional %s.\n"
			"   HINT: manually override with a scale factor, if known.\n", exCorrName.c_str());
	return iter->second;
}

VanDerWaals::AtomParams::AtomParams(double SI_C6, double SI_R0)
: C6(SI_C6 * Joule*pow(1e-9*meter,6)/mol), R0(SI_R0 * Angstrom)
{
}

VanDerWaals::AtomParams VanDerWaals::getParams(int atomicNumber) const
{	assert(atomicNumber>0);
	assert(atomicNumber<=atomicNumberMax);
	return atomParams[atomicNumber];
}


VanDerWaals::~VanDerWaals()
{
	for(auto& iter: radialFunctions) iter.second.free(); //Cleanup cached RadialFunctionG's if any
}

bool hackedTS = false;

struct HackedTSparams
{	double C6;
	double alpha;
	double R0;
	
	HackedTSparams(double C6, double C6free, double alphaFree, double R0free) : C6(C6)
	{	double Vratio = sqrt(C6/C6free);
		alpha = alphaFree * Vratio;
		R0 = R0free * pow(Vratio, 1./3);
	}

	static HackedTSparams get(int atomicNumber)
	{	switch(atomicNumber)
		{	case 1: return HackedTSparams(2.1, 6.5, 4.5, 1.64*Angstrom);
			case 6: return HackedTSparams(24.1, 46.6, 12, 1.90*Angstrom);
			case 7: return HackedTSparams(17.1, 24.2, 7.4, 1.77*Angstrom);
			case 8: return HackedTSparams(11.7, 15.6, 5.4, 1.66*Angstrom);
			case 16: return HackedTSparams(113, 134, 19.6, 2.04*Angstrom);
			case 17: return HackedTSparams(88.8, 94.6, 15, 1.96*Angstrom);
			default:
				die("No hacked TS params available.\n");
				return HackedTSparams(0, 0, 0, 0);
		}
	}
};

const RadialFunctionG& VanDerWaals::getRadialFunction(int atomicNumber1, int atomicNumber2) const
{
	//Check for precomputed radial function:
	std::pair<int,int> atomicNumberPair(
		std::min(atomicNumber1, atomicNumber2),
		std::max(atomicNumber1, atomicNumber2)); //Optimize Z1 <-> Z2 symmetry
	auto radialFunctionIter = radialFunctions.find(atomicNumberPair);
	if(radialFunctionIter != radialFunctions.end() && !hackedTS)
		return radialFunctionIter->second;
	
	//Get parameters for current pair of species:
	const AtomParams& params1 = getParams(atomicNumber1);
	const AtomParams& params2 = getParams(atomicNumber2);
	double C6 = sqrt(params1.C6 * params2.C6);
	double R0 = params1.R0 + params2.R0;
	
	if(hackedTS)
	{	/*const HackedTSparams& params1 = HackedTSparams::get(atomicNumber1);
		const HackedTSparams& params2 = HackedTSparams::get(atomicNumber2);
		C6 = 2.* params1.C6 * params2.C6 / (params1.C6*params2.alpha/params1.alpha + params2.C6*params1.alpha/params2.alpha);
		R0 = 0.1*(params1.R0 + params2.R0);*/
		R0 *= 0.1;
	}
	
	//Initialize function on real-space logarithmic radial grid:
	const double rMin = 1e-2;
	const double rMax = 1e+3;
	const double dlogr = 0.01;
	size_t nSamples = ceil(log(rMax/rMin)/dlogr);
	RadialFunctionR func(nSamples);
	double r = rMin, rRatio = exp(dlogr), E_r;
	for(size_t i=0; i<nSamples; i++)
	{	func.r[i] = r; //radial position
		func.dr[i] = r * dlogr; //integration weight
		func.f[i] = vdwPairEnergyAndGrad(r, C6, R0, E_r); //sample value
		r *= rRatio;
	}
	
	//Transform to reciprocal space, cache and return:
	RadialFunctionG& funcTilde = ((VanDerWaals*)this)->radialFunctions[atomicNumberPair];
	const double dGloc = 0.02; //same as the default for SpeciesInfo
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dGloc))+5;
	func.transform(0, dGloc, nGridLoc, funcTilde);
	return funcTilde;
}
