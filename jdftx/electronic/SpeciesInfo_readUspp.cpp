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

//! Length in bytes of Record markers used in Fortran sequential binary files
//! Technically depends on the Fortran compiler, but it is 4 for gfortran as well as ifort,
//! which covers what (almost) everyone uses to compile the USPP generator program.
#define FORTRAN_RECORD_MARKER_LENGTH 4 


//! Read the fortran sequential binary uspp format (need to skip record start/stop markers)
struct UsppReader
{	istream& is;
	char* recBuf;
	
	UsppReader(istream& is)
	: is(is)
	{	recBuf = new char[FORTRAN_RECORD_MARKER_LENGTH];
	}
	~UsppReader()
	{	delete[] recBuf;
	}
	
	//!Fortran record start
	void newRecord()
	{	is.read(recBuf, FORTRAN_RECORD_MARKER_LENGTH);  if(is.eof()) die("  File ended prematurely.\n");
	}
	//!Fortran record end
	void endRecord()
	{	is.read(recBuf, FORTRAN_RECORD_MARKER_LENGTH);  if(is.eof()) die("  Error reading file.\n");
	}
	
	//!Get a single element of type T
	template<typename T> void get(T& t)
	{	is.read((char*)&t, sizeof(T));
		if(is.eof()) die("  Error reading file.\n");
	}
	//!Get an array of length N of type T
	template<typename T> void get(T* t, size_t N)
	{	is.read((char*)t, sizeof(T)*N);
		if(is.eof()) die("  Error reading file.\n");
	}
	//!Get an vector of length N of type T
	template<typename T> void get(std::vector<T>& tVec)
	{	get(&tVec[0], tVec.size());
	}
	//!Get a string of length N
	template<int N> void get(string& s)
	{	char buf[N]; get(buf, N);
		s.assign(buf, N);
	}
};


void SpeciesInfo::readUspp(istream& is)
{
	UsppReader reader(is);
	
	//Header:
	reader.newRecord();
	int version[3]; reader.get(version, 3);
	int date[3]; reader.get(date, 3);
	reader.endRecord();
	
	reader.newRecord();
	string title; reader.get<20>(title); trim(title);
	double Zae; reader.get(Zae); //total electron count
	atomicNumber = int(round(Zae));
	reader.get(Z); //valence electron count
	double excorrCode; reader.get(excorrCode); //exchange-correlation code
	int nValence; reader.get(nValence); //number of valence orbitals
	int nGrid; reader.get(nGrid); //number of grid points in log mesh
	double Etot; reader.get(Etot); Etot *= 0.5; //total energy of atom (convert Ryd to Eh)
	reader.endRecord();
	
	logPrintf("  Title: %s.  Created by USPP %d.%d.%d on %d-%d-%d\n",
		title.c_str(), version[0], version[1], version[2], date[0], date[1], date[2]);
	logPrintf("  Reference state energy: %lf.  %lg valence electrons in orbitals:\n", Etot, Z);
	if(version[0]<7 || (version[0]==7 && version[1]<3))
		die("  Only USPP format 7.2 or newer is supported.\n");
	
	//Valence orbital (vo) properties
	reader.newRecord();
	std::vector<int> voCode(nValence); //n,l,m
	std::vector<double> voWeight(nValence); //fillings in ref state
	std::vector<double> voEnergy(nValence); //eigenvalue
	for(int v=0; v<nValence; v++)
	{	reader.get(voCode[v]);
		reader.get(voWeight[v]);
		reader.get(voEnergy[v]); voEnergy[v] *= 0.5; //convert Ryd to Eh
		logPrintf("    |%d>  occupation: %lg  eigenvalue: %lf\n", voCode[v], voWeight[v], voEnergy[v]);
	}
	reader.endRecord();
	
	//pseudopotential generation method (ignored), core presence
	reader.newRecord();
	int keyPS; reader.get(keyPS);
	int haveCore; reader.get(haveCore); //non-zero if core correction present
	double rInnerIgnored; reader.get(rInnerIgnored);
	reader.endRecord();
	
	//Number of angular channels, Qij coefficients etc:
	reader.newRecord();
	int lMax; reader.get(lMax); lMax--;  //maximum projetcor angular momentum
	int lLocal; reader.get(lLocal); //local channel (<0 if custom)
	double Elocal; reader.get(Elocal); Elocal *= 0.5; //local energy
	int ifqOpt; reader.get(ifqOpt); //some flag to do with Qij pseudization (irrelevant here)
	int nCoeff; reader.get(nCoeff); //number of Taylor series coefficients for Qij within rInner
	double QijEcut; reader.get(QijEcut); QijEcut *= 0.5; //energy cutoff used for Qij pseudization
	reader.endRecord();
	logPrintf("  lMax: %d  lLocal: %d  QijEcut: %lg\n", lMax, lLocal, QijEcut);
	if(lMax>3) die("  Nonlocal projectors with l>3 not implemented (lMax = %d not supported).\n", lMax);
	
	int nL = 2*lMax+1; //number of l-channels in pair products of Ylm's (0 to 2 lMax)
	reader.newRecord();
	std::vector<double> rInner(nL); reader.get(rInner); //radius upto which Qij is pseudized
	reader.endRecord();
	
	reader.newRecord();
	int iRel; reader.get(iRel); //relativity type used in pseudopotential generation
	reader.endRecord();

	reader.newRecord();
	std::vector<double> rcNL(lMax+1); reader.get(rcNL); //non-local projector radius
	reader.endRecord();

	reader.newRecord();
	int nBeta; reader.get(nBeta);
	int nGridBeta; reader.get(nGridBeta);
	reader.endRecord();
	logPrintf("  %d projectors sampled on a log grid with %d points:\n", nBeta, nGridBeta);
	
	std::vector<int> lNL(nBeta); //l for each projector
	std::vector<std::vector<int> > lBeta(lMax+1); //Projector indices for each l
	std::vector<double> eigNL(nBeta); //eigenvalue for each projector
	std::vector<RadialFunctionR> Vnl(nBeta, RadialFunctionR(nGridBeta)); //projector radial function
	std::vector<std::vector<double> > D0(nBeta,std::vector<double>(nBeta)); //descreened D matrix (Mnl in our notation)
	std::vector<std::vector<double> > D(nBeta,std::vector<double>(nBeta)); //full D matrix (not needed)
	std::vector<std::vector<double> > Q(nBeta,std::vector<double>(nBeta)); //Q matrix (goes to Qint in SpeciesInfo)
	std::vector<std::vector<RadialFunctionR> > Qr(nBeta, std::vector<RadialFunctionR>(nBeta, RadialFunctionR(nGridBeta))); //overlap radial function
	std::vector<std::vector<std::vector<std::vector<double> > > > Qcoeff(nBeta, std::vector<std::vector<std::vector<double> > >(nBeta));
	for(int i=0; i<nBeta; i++)
	{	reader.newRecord();
		reader.get(lNL[i]); lBeta[lNL[i]].push_back(i); //Update index map
		reader.get(eigNL[i]); eigNL[i] *= 0.5; //Ryd to Eh
		reader.get(Vnl[i].f);
		reader.endRecord();
		logPrintf("    l: %d  eig: %lf  rCut: %lg\n", lNL[i], eigNL[i], rcNL[lNL[i]]);
		
		for(int j=i; j<nBeta; j++)
		{	Qcoeff[i][j].resize(nL, std::vector<double>(nCoeff));
			reader.newRecord();
			reader.get(D0[i][j]);
			reader.get(D[i][j]);
			reader.get(Q[i][j]);
			reader.get(Qr[i][j].f);
			for(int iL=0; iL<nL; iL++)
				reader.get(Qcoeff[i][j][iL]);
			reader.endRecord();
		}
	}

	reader.newRecord();
	std::vector<int> ipType(nBeta); reader.get(ipType);
	int nP; reader.get(nP);
	double ptryc; reader.get(ptryc);
	reader.endRecord();

	//Descreened local potential:
	reader.newRecord();
	double rcLoc; reader.get(rcLoc);
	RadialFunctionR Vloc0(nGrid); reader.get(Vloc0.f);
	reader.endRecord();

	double rcCore=0.; RadialFunctionR nCore(nGrid);
	if(haveCore)
	{	reader.newRecord(); reader.get(rcCore); reader.endRecord();
		reader.newRecord(); reader.get(nCore.f); reader.endRecord();
		logPrintf("  Partial core density with radius %lg\n", rcCore);
	}

	//Local potential in reference configuration (not used)
	reader.newRecord();
	RadialFunctionR Vloc(nGrid); reader.get(Vloc.f);
	reader.endRecord();
	
	//Charge density
	reader.newRecord();
	RadialFunctionR rsAtom(nGrid); reader.get(rsAtom.f);
	reader.endRecord();

	//Radial locations:
	reader.newRecord();
	std::vector<double> rGrid(nGrid); reader.get(rGrid);
	reader.endRecord();

	//Radial weights:
	reader.newRecord();
	std::vector<double> drGrid(nGrid); reader.get(drGrid);
	reader.endRecord();

	//Wavefunctions:
	reader.newRecord();
	int nPsi; reader.get(nPsi);
	reader.endRecord();
	std::vector<RadialFunctionR> voPsi(nPsi, RadialFunctionR(nGrid)); //radial wavefunctions
	reader.newRecord();
	for(int v=0; v<nPsi; v++) reader.get(voPsi[v].f);
	reader.endRecord();
	
	//-------------- Transform and store required quantities in SpeciesInfo ------------
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	int nGridNL = int(ceil(e->gInfo.GmaxSphere/dG))+5;
	
	//Core density:
	if(haveCore)
	{	nCore.set(rGrid, drGrid);
		for(int i=0; i<nGrid; i++)
			nCore.f[i] *= (rGrid[i] ? 1./(4*M_PI*rGrid[i]*rGrid[i]) : 0);
		setCore(nCore);
	}
	
	//Local potential:
	Vloc0.set(rGrid, drGrid);
	logPrintf("  Transforming local potential to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridLoc);
	for(int i=0; i<nGrid; i++)
		Vloc0.f[i] = (Vloc0.f[i]*0.5 + Z) * (rGrid[i] ? 1./rGrid[i] : 0); //Convert to Eh and remove the -Z/r part
	Vloc0.transform(0, dG, nGridLoc, VlocRadial);
	
	//Projectors:
	if(nBeta)
	{	logPrintf("  Transforming nonlocal projectors to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridNL);
		VnlRadial.resize(lMax+1);
		for(int iBeta=0; iBeta<nBeta; iBeta++)
		{	int l = lNL[iBeta];
			VnlRadial[l].push_back(RadialFunctionG());
			Vnl[iBeta].set(rGrid, drGrid);
			for(int i=0; i<nGridBeta; i++)
				Vnl[iBeta].f[i] *= (rGrid[i] ? 1./rGrid[i] : 0);
			Vnl[iBeta].transform(l, dG, nGridNL, VnlRadial[l].back());
		}
		//Set Mnl:
		Mnl.resize(lMax+1);
		for(int l=0; l<=lMax; l++)
		{	int nProj = lBeta[l].size();
			Mnl[l].init(nProj, nProj);
			for(int p1=0; p1<nProj; p1++)
				for(int p2=0; p2<nProj; p2++)
				{	int iBeta = lBeta[l][p1];
					int jBeta = lBeta[l][p2];
					if(iBeta>jBeta) std::swap(iBeta,jBeta);
					Mnl[l].data()[Mnl[l].index(p1,p2)] = D0[iBeta][jBeta]*0.5; //Convert to Eh
				}
		}
		//Q radial functions and integrals:
		logPrintf("  Transforming density augmentations to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridLoc);
		Qint.resize(lMax+1);
		for(int l=0; l<=lMax; l++) if(lBeta[l].size())
			Qint[l] = zeroes(lBeta[l].size(), lBeta[l].size());
		for(int l1=0; l1<=lMax; l1++) for(int p1=0; p1<int(lBeta[l1].size()); p1++)
		{	int iBeta = lBeta[l1][p1];
			for(int l2=0; l2<=lMax; l2++) for(int p2=0; p2<int(lBeta[l2].size()); p2++)
			{	int jBeta = lBeta[l2][p2];
				if(iBeta<=jBeta) //store and use upper triangular only
				{	RadialFunctionR& Qij = Qr[iBeta][jBeta];
					for(int l=abs(l1-l2); l<=l1+l2; l+=2)
					{	//Replace the inner section with the l-dependent pseudization:
						std::vector<double>& coeff = Qcoeff[iBeta][jBeta][l];
						RadialFunctionR Qijl(nGridBeta); Qijl.set(rGrid, drGrid);
						bool isNonzero = false; int iLastNZ=0; double Qabs = 0.;
						for(int i=0; i<nGridBeta; i++) 
						{	if(rGrid[i]<rInner[l])
							{	double val=0.0, rSq=rGrid[i]*rGrid[i], rPow=pow(rGrid[i],l);
								for(auto c: coeff)
								{	val += c*rPow;
									rPow *= rSq;
								}
								Qijl.f[i] = val;
							}
							else Qijl.f[i] = Qij.f[i]/pow(rGrid[i],2);
							//Radial extent check:
							double dQabs = fabs(Qijl.f[i]) * rGrid[i]*rGrid[i] * drGrid[i];
							Qabs += dQabs;
							if(dQabs>1e-16*Qabs) { isNonzero=true; iLastNZ=i+2; }
						}
						if(!isNonzero) continue;
						if(iLastNZ < nGridBeta)
						{	Qijl.r.resize(iLastNZ);
							Qijl.dr.resize(iLastNZ);
							Qijl.f.resize(iLastNZ);
						}
						//Store in Qradial:
						QijIndex qIndex = { l1, p1, l2, p2, l };
						Qijl.transform(l, dG, nGridLoc, Qradial[qIndex]);
						//Store Qint = integral(Qradial) when relevant:
						if(l1==l2 && !l)
						{	double Qint_ij = Qijl.transform(0,0)/(4*M_PI);
							Qint[l1].data()[Qint[l1].index(p1,p2)] = Qint_ij;
							Qint[l1].data()[Qint[l1].index(p2,p1)] = Qint_ij;
						}
					}
				}
			}
		}
		if(!Qradial.size()) Qint.clear(); //Special case of norm-conserving PSP in USPP skin
	}
	
	//Wavefunctions:
	if(nPsi == nValence) // <- this seems to always be the case, but just to be sure
	{	logPrintf("  Transforming atomic orbitals to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridNL);
		psiRadial.resize(lMax+1);
		atomEigs.resize(lMax+1);
		for(int v=0; v<nValence; v++)
		{	int l = (voCode[v]/10)%10; //l is 2nd digit of orbital code (ie. |210> is l=1)
			if(l>lMax) //this can happen occassionaly when lLocal==lMax
			{	psiRadial.resize(l+1);
				atomEigs.resize(l+1);
			}
			voPsi[v].set(rGrid, drGrid);
			for(int i=0; i<nGrid; i++)
				voPsi[v].f[i] *= (rGrid[i] ? 1./rGrid[i] : 0);
			psiRadial[l].push_back(RadialFunctionG());
			voPsi[v].transform(l, dG, nGridNL, psiRadial[l].back());
			//Store eigenvalue:
			atomEigs[l].push_back(voEnergy[v]);
		}
	}
	
	//Determine max core radius:
	std::set<double> coreRadii; //collection of all relevant radii that should not be allowed to overlap
	coreRadii.insert(rInner.begin(), rInner.end()); //Qij pesudization radius
	coreRadii.insert(rcNL.begin(), rcNL.end()); //Non-local projector radius
	coreRadii.insert(rcLoc); //Local channel pseudization radius
	coreRadii.insert(rcCore); //Partial-core pseudization radius
	coreRadius = *coreRadii.rbegin(); //max of all radii above
}
