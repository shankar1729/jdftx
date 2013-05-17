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

#include <electronic/SpeciesInfo.h>
#include <electronic/Everything.h>
#include <electronic/matrix.h>


namespace PotFile
{	
	//! Return the first non-blank, non-comment line:
	string getLine(istream& in)
	{	string line;
		while(line.find_first_not_of(" \t\n\r")==string::npos || line.at(0)=='#')
			getline(in, line);
		return line;
	}

	//! Decsription of a radial function
	struct RadialInfo
	{	int nPoints;
		double dG;
		string fname;
		
		void readInfo(istream& in, bool local, double Gmax)
		{	istringstream iss(PotFile::getLine(in));
			iss >> nPoints >> dG >> fname;
			logPrintf("nPoints: %d   dG: %g  file: %s\n", nPoints, dG, fname.c_str());
			//Check if q-grid is sufficient:
			if((nPoints-5)*dG < Gmax)
			{	logPrintf("        WARNING: Insufficient resolution for this %s:\n", local ? "FFT box" : "Ecut");
				logPrintf("        (nPoints-5)*dG = %lg < Gmax%s = %lg\n", (nPoints-5)*dG, local?"Loc":"NL", Gmax);
			}
		}
		
		void readFunction(int l, RadialFunctionG& func, double scale=1.0) const
		{	func.init(l, nPoints, dG, fname.c_str(), scale);
		}
		
		bool operator==(const RadialInfo& other) const
		{	return (nPoints==other.nPoints) && (dG==other.dG) && (fname==other.fname);
		}
	};

	//! Information for a given angular momentum state (l,m)
	struct AngularInfo
	{	int l; //!< angular quantum numbers
		std::vector<RadialInfo> proj; //!< projector info
		matrix M; //!< projector matrix
		
		void readInfoVnl(istream& in, double GmaxNL)
		{	int m; //read in but not used
			int nProj;
			istringstream(PotFile::getLine(in)) >> l >> m >> nProj;
			logPrintf("    l: %d  m: %d  nProj: %d\n", l, m, nProj);
			//Matrix:
			M.init(nProj, nProj);
			complex* Mdata = M.data();
			for(int i=0; i<nProj; i++)
				for(int j=0; j<nProj; j++)
				{	complex& c = Mdata[M.index(i,j)];
					in >> c.real();
					c.imag() = 0;
				}
			logPrintf("      M: \n");
			M.print_real(globalLog, "        %lg");
			//Radial functions:
			logPrintf("      Radial functions:\n");
			proj.resize(nProj);
			for(RadialInfo& p: proj) { logPrintf("        "); p.readInfo(in, false, GmaxNL); }
		}
		
		void readInfoPAW(istream& in, double GmaxNL)
		{	int nProj;
			istringstream(PotFile::getLine(in)) >> l >> nProj;
			logPrintf("    l: %d  nProj: %d\n", l, nProj);
			logPrintf("      Radial functions:\n");
			proj.resize(nProj);
			for(RadialInfo& p: proj) { logPrintf("        "); p.readInfo(in, false, GmaxNL); }
		}
		
		void readFunctions(std::vector<RadialFunctionG>& funcs) const
		{	funcs.resize(proj.size());
			for(unsigned p=0; p<proj.size(); p++)
				proj[p].readFunction(l, funcs[p], 4*M_PI); // Internal normalization is different from legacy POT file format
		}
		
		bool operator!=(const AngularInfo& other)
		{	return (l!=other.l) || (proj!=other.proj) || nrm2(M-other.M);
		}
	};
}

void SpeciesInfo::readPot(istream& in)
{
	using namespace PotFile;

	//Local potential:
	logPrintf("  Local potential:  ");
	RadialInfo VlocInfo;
	VlocInfo.readInfo(in, true, e->gInfo.GmaxGrid);
	VlocInfo.readFunction(0, VlocRadial);
	
	//Non-local potential:
	int nlm; istringstream(getLine(in)) >> nlm;
	logPrintf("  Non-local potential:  nlm: %d\n", nlm);
	if(nlm != 0)
	{	//First read by combined lm index:
		std::map<int, AngularInfo> lInfoMap; //map from angular momentum l to info
		int lMax = 0;
		for(int lm=0; lm<nlm; lm++)
		{	AngularInfo temp; temp.readInfoVnl(in, e->gInfo.GmaxSphere);
			if(temp.l > lMax) lMax = temp.l;
			if(temp.l < 0) die("  Encountered l<0 in pseudopotential.\n");
			if(lInfoMap.find(temp.l)==lInfoMap.end()) lInfoMap[temp.l] = temp; //new l
			else //existing l, make sure info agrees
				if(temp != lInfoMap[temp.l])
					die("  Pseudopotential breaks spherical symmetry at l=%d (info is m-dependent).\n", temp.l);
		}
		if(lMax>3) die("  Nonlocal projectors with l>3 not implemented (lMax = %d not supported).\n", lMax);
		//Now read the radial functions for each l (have ensured no m dependence)
		VnlRadial.resize(lMax+1);
		Mnl.resize(lMax+1);
		for(auto lInfoEntry: lInfoMap)
		{	const int& l = lInfoEntry.first;
			const AngularInfo& lInfo = lInfoEntry.second;
			lInfo.readFunctions(VnlRadial[l]);
			Mnl[l] = lInfo.M;
		}
	}
}

// Read pulay stuff from file which has a line with number of Ecuts
// followed by arbitrary number of Ecut dE_dnG pairs
void SpeciesInfo::setupPulay()
{	
	using namespace PotFile;
	if(pulayfilename == "none") return;
	
	ifstream ifs(pulayfilename.c_str());
	if(!ifs.is_open()) die("  Can't open pulay file %s for reading.\n", pulayfilename.c_str());
	logPrintf("  Reading pulay file %s ... ", pulayfilename.c_str());
	istringstream iss;
	int nEcuts; istringstream(getLine(ifs)) >> nEcuts;
	std::map<double,double> pulayMap;
	for(int i=0; i<nEcuts; i++)
	{	double Ecut, dE_dnG;
		istringstream(getLine(ifs)) >> Ecut >> dE_dnG;
		pulayMap[Ecut] = dE_dnG;
	}
	
	double minEcut = pulayMap.begin()->first;
	double maxEcut = pulayMap.rbegin()->first;
	double Ecut = e->cntrl.Ecut;
	if(pulayMap.find(Ecut) != pulayMap.end())
	{	dE_dnG = pulayMap[Ecut];
		logPrintf("using dE_dnG = %le computed for Ecut = %lg.\n", dE_dnG, Ecut);
	}
	else if(Ecut < minEcut)
	{	die("\n  Ecut=%lg < smallest Ecut=%lg in pulay file %s.\n",
			Ecut, minEcut, pulayfilename.c_str());
	}
	else if(Ecut > maxEcut)
	{	dE_dnG = 0;
		logPrintf("using dE_dnG = 0 for Ecut=%lg > largest Ecut=%lg in file.\n",
			Ecut, maxEcut);
	}
	else
	{	auto iRight = pulayMap.upper_bound(Ecut); //iterator just > Ecut
		auto iLeft = iRight; iLeft--; //iterator just < Ecut
		double t = (Ecut - iLeft->first) / (iRight->first - iLeft->first);
		dE_dnG = (1.-t) * iLeft->second + t * iRight->second;
		logPrintf("using dE_dnG = %le interpolated from Ecut = %lg and %lg.\n",
			dE_dnG, iLeft->first, iRight->first);
	}
}
