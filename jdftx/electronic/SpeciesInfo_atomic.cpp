/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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
#include <electronic/SpeciesInfo_internal.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/LatticeUtils.h>

//------- SpeciesInfo functions related to atomic orbitals -------


//Overlap of two radial functions (assumes same grid configuration, but one grid could be shorter)
double dot(const RadialFunctionR& X, const RadialFunctionR& Y)
{	size_t nr = std::min(X.f.size(), Y.f.size());
	assert(X.r.size() >= nr);
	assert(X.dr.size() >= nr);
	double ret = 0.;
	for(size_t i=0; i<nr; i++)
	{	const double& r = X.r[i];
		const double& dr = X.dr[i];
		ret += (r*r*dr) * (X.f[i] * Y.f[i]);
	}
	return ret;
}
//Accumulate radial functions (assumes same grid configuration, but X could be shorter)
void axpy(double alpha, const RadialFunctionR& X, RadialFunctionR& Y)
{	size_t nr = X.f.size();
	assert(Y.f.size() >= nr);
	for(size_t i=0; i<nr; i++) Y.f[i] += alpha * X.f[i];
}



void SpeciesInfo::accumulateAtomicDensity(ScalarFieldTildeArray& nTilde) const
{	//Collect list of distinct magnetizations and corresponding atom positions:
	struct MomentDirection //list of atoms with same magnetic moment
	{	vector3<> Mhat;
		std::vector< vector3<> > atpos;
	};
	struct MomentMagnitude //list of atoms with same magnetic moment magnitude
	{	double Mlength;
		std::vector<MomentDirection> Mdirs;
	};
	std::vector<MomentMagnitude>  Mmags;
	if(initialMagneticMoments.size())
	{	for(unsigned a=0; a<atpos.size(); a++)
		{	const vector3<>& M = initialMagneticMoments[a];
			double Mlength = M.length();
			vector3<> Mhat = (Mlength ? 1./Mlength : 0.) * M;
			bool foundLength = false;
			for(MomentMagnitude& Mmag: Mmags) if(fabs(Mlength - Mmag.Mlength) < symmThreshold)
			{	bool foundDir = false;
				for(MomentDirection& Mdir: Mmag.Mdirs) if((M - Mmag.Mlength*Mdir.Mhat).length_squared() < symmThresholdSq)
				{	Mdir.atpos.push_back(atpos[a]);
					foundDir = true;
					break;
				}
				if(!foundDir)
				{	MomentDirection Mdir;
					Mdir.Mhat = Mhat;
					Mdir.atpos.push_back(atpos[a]);
					Mmag.Mdirs.push_back(Mdir);
				}
				foundLength = true;
				break;
			}
			if(!foundLength)
			{	MomentDirection Mdir;
				Mdir.Mhat = Mhat;
				Mdir.atpos.push_back(atpos[a]);
				MomentMagnitude Mmag;
				Mmag.Mlength = Mlength;
				Mmag.Mdirs.push_back(Mdir);
				Mmags.push_back(Mmag);
			}
		}
	}
	else
	{	MomentDirection Mdir;
		Mdir.atpos = atpos;
		MomentMagnitude Mmag;
		Mmag.Mlength = 0.;
		Mmag.Mdirs.push_back(Mdir);
		Mmags.push_back(Mmag);
	}
	//For each magnetization magnitude:
	for(const MomentMagnitude& Mmag: Mmags)
	{	if(e->eInfo.nDensities == 1) assert(!Mmag.Mlength);
		//Get the radial spin-densities in diagonal basis:
		std::vector<RadialFunctionG> nRadial(Mmag.Mlength ? 2 : 1);
		for(unsigned s=0; s<nRadial.size(); s++)
			getAtom_nRadial(s, Mmag.Mlength, nRadial[s], false);
		//Collect contributions for each direction with this magnitude:
		for(const MomentDirection& Mdir: Mmag.Mdirs)
		{	//Compute structure factor for atoms with current magnetization vector:
			ManagedArray<vector3<>> atposCur(Mdir.atpos);
			ScalarFieldTilde SG; nullToZero(SG, e->gInfo);
			callPref(getSG)(e->gInfo.S, Mdir.atpos.size(), atposCur.dataPref(), 1./e->gInfo.detR, SG->dataPref());
			//Spin-densities in diagonal basis
			std::vector<ScalarFieldTilde> nDiag;
			for(const RadialFunctionG& nRad: nRadial)
				nDiag.push_back(nRad * SG);
			//Accumulate contributions:
			if(nDiag.size()==1)
			{	if(nTilde.size()==1)
					nTilde[0] += 2.*nDiag[0];
				else
				{	nTilde[0] += nDiag[0];
					nTilde[1] += nDiag[0];
				}
			}
			else
			{	nTilde[0] += (0.5*(1.+Mdir.Mhat[2])) * nDiag[0] + (0.5*(1.-Mdir.Mhat[2])) * nDiag[1];
				nTilde[1] += (0.5*(1.-Mdir.Mhat[2])) * nDiag[0] + (0.5*(1.+Mdir.Mhat[2])) * nDiag[1];
				if(nTilde.size()==4) //noncollinear
				{	nTilde[2] += (0.5*Mdir.Mhat[0]) * (nDiag[0] - nDiag[1]); //Re(UpDn) = 0.5*Mx
					nTilde[3] += (0.5*Mdir.Mhat[1]) * (nDiag[1] - nDiag[0]); //Im(UpDn) = -0.5*My
				}
			}
		}
		for(RadialFunctionG& nRad: nRadial) nRad.free();
	}
}

void SpeciesInfo::accumulateAtomicPotential(ScalarFieldTilde& dTilde) const
{	//Get radial potential due to one atom:
	RadialFunctionG dRadial;
	getAtomPotential(dRadial);
	//Gets tructure factor:
	ScalarFieldTilde SG; nullToZero(SG, e->gInfo);
	callPref(getSG)(e->gInfo.S, atpos.size(), atposManaged.dataPref(), 1./e->gInfo.detR, SG->dataPref());
	//Accumulate contrbutions:
	dTilde += dRadial * SG;
	dRadial.free();
}

//Set atomic orbitals in column bundle from radial functions (almost same operation as setting Vnl)
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& Y, bool applyO, int colOffset, const vector3<>* derivDir) const
{	if(!atpos.size()) return;
	const auto& fRadial = applyO ? OpsiRadial : psiRadial; //!< select radial function set (psi or Opsi)
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	int nOrbitalsPerAtom = 0;
	for(int l=0; l<int(fRadial.size()); l++)
		nOrbitalsPerAtom += nAtomicOrbitals(l)*(2*l+1)*nSpinCopies;
	int iCol = colOffset;
	for(int l=0; l<int(fRadial.size()); l++)
		for(int n=0; n<nAtomicOrbitals(l); n++)
		{	setAtomicOrbitals(Y, applyO, n, l, iCol, nOrbitalsPerAtom, derivDir);
			iCol += (2*l+1)*nSpinCopies;
		}
}
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& psi, bool applyO, unsigned n, int l, int colOffset, int atomColStride, const vector3<>* derivDir) const
{	if(!atpos.size()) return;
	assert(l < int(psiRadial.size()));
	assert(int(n) < nAtomicOrbitals(l));
	const auto& fRadial = applyO ? OpsiRadial : psiRadial; //!< select radial function set (psi or Opsi)
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	int nOrbitalsPerAtom = (2*l+1)*nSpinCopies;
	if(atomColStride) assert(atomColStride >= nOrbitalsPerAtom); else atomColStride = nOrbitalsPerAtom;
	assert(psi.basis); assert(psi.qnum);
	assert(colOffset + atomColStride*int(atpos.size()-1) + nOrbitalsPerAtom <= psi.nCols());
	if(nSpinCopies>1) assert(psi.isSpinor()); //can have multiple spinor copies only in spinor mode
	const Basis& basis = *psi.basis;
	if(isRelativistic() && l>0)
	{	//find the two orbital indices corresponding to different j of same n
		std::vector<int> pArr; 
		for(int j2=2*l-1; j2<=2*l+1; j2+=2)
		{	unsigned ns = 0;
			for(unsigned p=0; p<psiRadial[l].size(); p++)
				if(psi2j[l][p]==j2)
				{	if(ns==n) { pArr.push_back(p); break; }
					ns++;
				}
		}
		//Initialize a non-spinor ColumnBundle containing all m's for both j functions:
		ColumnBundle V(atpos.size()*nOrbitalsPerAtom, basis.nbasis, &basis, psi.qnum, isGpuEnabled());
		int iCol=0;
		for(int p: pArr) for(int m=-l; m<=l; m++)
		{	size_t atomStride = V.colLength() * nOrbitalsPerAtom;
			size_t offs = iCol * V.colLength();
			callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, psi.qnum->k, basis.iGarr.dataPref(),
				e->gInfo.G, atposManaged.dataPref(), fRadial[l][p], V.dataPref()+offs, derivDir);
			iCol++;
		}
		//Transform the non-spinor ColumnBundle to the spinorial j eigenfunctions:
		int N = nOrbitalsPerAtom;
		matrix transform = zeroes(2*N, N);
		transform.set(0,N,   0,2*l, getYlmToSpinAngleMatrix(l, 2*l-1));
		transform.set(N,2*N, 2*l,N, getYlmToSpinAngleMatrix(l, 2*l+1));
		for(size_t a=0; a<atpos.size(); a++)
			psi.setSub(colOffset+a*atomColStride, V.getSub(a*N,(a+1)*N) * transform);
	}
	else
	{	int iCol = colOffset; //current column
		for(int m=-l; m<=l; m++)
		{	//Set atomic orbitals for all atoms at specified (n,l,m):
			size_t atomStride = psi.colLength() * atomColStride;
			size_t offs = iCol * psi.colLength();
			callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, psi.qnum->k, basis.iGarr.dataPref(),
				e->gInfo.G, atposManaged.dataPref(), fRadial[l][n], psi.dataPref()+offs, derivDir);
			if(nSpinCopies>1) //make copy for other spin
			{	complex* dataPtr = psi.dataPref()+offs;
				for(size_t a=0; a<atpos.size(); a++)
				{	callPref(eblas_zero)(2*basis.nbasis, dataPtr+basis.nbasis);
					callPref(eblas_copy)(dataPtr+3*basis.nbasis, dataPtr, basis.nbasis);
					dataPtr += atomStride;
				}
			}
			iCol += nSpinCopies;
		}
	}
}
int SpeciesInfo::nAtomicOrbitals() const
{	int nOrbitals = 0;
	if(isRelativistic())
	{	for(int l=0; l<int(psiRadial.size()); l++)
			for(unsigned n=0; n<psiRadial[l].size(); n++)
				nOrbitals += (psi2j[l][n]+1);
		return nOrbitals * atpos.size();
	}
	else
	{	for(int l=0; l<int(psiRadial.size()); l++)
			nOrbitals += (2*l+1)*psiRadial[l].size();
		int nSpinCopies = 2/e->eInfo.qWeightSum;
		return nOrbitals * atpos.size() * nSpinCopies;
	}
}
int SpeciesInfo::lMaxAtomicOrbitals() const
{	return int(psiRadial.size()) - 1;
}
int SpeciesInfo::nAtomicOrbitals(int l) const
{	assert(l >= 0);
	if(unsigned(l) >= psiRadial.size()) return -1; //signals end of l
	return psiRadial[l].size() / ((isRelativistic() && l>0) ? 2 : 1); //principal quantum number reduced by a factor of 2 when psiRadial counts the j splitting
}
int SpeciesInfo::atomicOrbitalOffset(unsigned int iAtom, unsigned int n, int l, int m, int s) const
{	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(s < e->eInfo.spinorLength());
	assert(int(n) < nAtomicOrbitals(l));
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	int iOrb = 0;
	for(int L=0; L<int(psiRadial.size()); L++)
	{	int nCount_L = nAtomicOrbitals(L);
		iOrb += nCount_L * (2*L+1)*nSpinCopies * iAtom; //orbitals from previous atoms
		if(L<=l)
			iOrb += (L==l ? n : nCount_L) * (2*L+1)*nSpinCopies; //orbitals from previous l,n at current atom
	}
	if(isRelativistic())
	{	int j2 = 2*l + (s ? -1 : +1); //2*j
		int mj2 = 2*m +  + (s ? -1 : +1); //2*m_j
		assert(mj2 >= -j2); assert(mj2 <= j2);
		if(s==0) iOrb += 2*l; //include orbitals from previous j at current n,l,atom (j = l-1/2 (accessed using s=1) is stored first)
		return iOrb + (j2+mj2)/2; //include orbitals from previous m at current j,n,l,atom
	}
	else
	{	assert(m >= -l); assert(m <= l);
		return iOrb + (l+m)*nSpinCopies + s; //include orbitals from previous m,s at current n,l,atom
	}
}

void SpeciesInfo::estimateAtomEigs()
{	if(!psiRadial.size()) return; //no orbitals

	//Compute eigenvalues if unavailable
	if(!atomEigs.size())
	{	std::map<int,int> invalidPsis; //some FHI pseudodpotentials store unbound projectors which need to be discarded for LCAO
		logPrintf("  Approximate pseudo-atom eigenvalues: ");
		atomEigs.resize(psiRadial.size());
		double normFac = 1./e->gInfo.detR; //normalization factor in wavefunctions
		for(unsigned l=0; l<psiRadial.size(); l++)
		{	for(unsigned p=0; p<psiRadial[l].size(); p++)
			{	const RadialFunctionR& psi = *(psiRadial[l][p].rFunc);
				//Find points deep in the tail of the wavefunction, but still far above roundoff limit:
				double r[2], e[2];
				double psiThresh[2] = { 3e-7*normFac, 3e-6*normFac };
				unsigned i = psi.r.size()-2;
				for(int j=0; j<2; j++)
				{	for(; i>1; i--) if(fabs(psi.f[i]) > psiThresh[j]) break;
					//Apply the Schrodinger operator with V = 0 (tail of a neutral atom)
					r[j] = psi.r[i]; double rm = psi.r[i-1], rp = psi.r[i+1];
					double Dfp = (psi.f[i+1]/psi.f[i] - 1.) / (rp-r[j]);
					double Dfm = (psi.f[i-1]/psi.f[i] - 1.) / (rm-r[j]);
					if(Dfp >= 0. || Dfm >= 0.) { e[j] = 0.; continue; } //must be an unbound state
					e[j] = (-0.5/(r[j]*r[j])) * (l*(l+1) + (0.5/(rp-rm)) * (Dfp*(rp+r[j])*(rp+r[j]) - Dfm*(rm+r[j])*(rm+r[j])));
				}
				if(e[0] && e[1])
				{	double eCorrected = (r[0]*e[0] - r[1]*e[1]) / (r[0] - r[1]); //correct for remnant Z/r term in non-neutral atoms
					atomEigs[l].push_back(eCorrected);
				}
				else //unbound state
				{	atomEigs[l].push_back(std::numeric_limits<double>::infinity()); //make sure it doesn't get occupied
					invalidPsis[l]++;
				}
			}
			//Print:
			const char orbCode[] = "spdfgh";
			if(isRelativistic() && l>0) //split by j
			{	int l2 = 2*l;
				for(int j2=l2-1; j2<=l2+1; j2+=2)
				{	logPrintf("  %c%c (", orbCode[l], (j2<l2 ? '-' : '+'));
					for(unsigned p=0; p<psiRadial[l].size(); p++)
						if(psi2j[l][p]==j2) logPrintf(" %.2lg", atomEigs[l][p]);
					logPrintf(" )");
				}
			}
			else
			{	logPrintf("  %c (", orbCode[l]);
				for(unsigned p=0; p<psiRadial[l].size(); p++)
					logPrintf(" %.2lg", atomEigs[l][p]);
				logPrintf(" )");
			}
		}
		logPrintf("\n");
		//Report removal of invalid psi's:
		for(auto invalidPsi: invalidPsis)
			logPrintf("  WARNING: encountered %d l=%d unbound projectors in atomic orbital set.\n", invalidPsi.second, invalidPsi.first);
	}
}

void SpeciesInfo::getAtom_nRadial(int spin, double magneticMoment, RadialFunctionG& nRadial, bool forceNeutral) const
{
	int spinCount = (e->eInfo.nDensities==1 ? 1 : 2);
	assert(spin >= 0); assert(spin < spinCount);
	
	//Determine occupations:
	std::multimap< double, std::pair<unsigned,unsigned> > eigMap; //map from eigenvalues to (l,p)
	std::vector<std::vector<double> > atomF(psiRadial.size());
	for(unsigned l=0; l<psiRadial.size(); l++)
	{	atomF[l].resize(psiRadial[l].size());
		for(unsigned p=0; p<psiRadial[l].size(); p++)
			eigMap.insert(std::make_pair(atomEigs[l][p], std::make_pair(l,p)));
	}
	if(spinCount>1 && magneticMoment)
		logPrintf("%s (M=%lg) pseudo-atom %s spin occupations: ", name.c_str(), magneticMoment, spin==0 ? "majority" : "minority");
	else
		logPrintf("%s pseudo-atom occupations: ", name.c_str());
	double Nvalence = Z - initialOxidationState;
	double N = 0.5*(Nvalence + (1-2*spin)*magneticMoment); //total electrons to fill
	if(forceNeutral) N = 0.5*Z; //force neutral unpolarized atom
	if(N < 0.)
		die("Magnetic moment (%lg) exceeds pseudo-atom valence electron count (%lg) [per spin channel].\n", magneticMoment, Nvalence);
	double Favail = N; //electrons yet to be filled
	for(auto eigEntry: eigMap) //in ascending order of eigenvalues
	{	unsigned l = eigEntry.second.first;
		unsigned p = eigEntry.second.second;
		double capacity = isRelativistic() ? 0.5*(psi2j[l][p]+1) : (2*l+1);
		atomF[l][p] = std::min(Favail, capacity);
		Favail -= atomF[l][p];
	}
	if(Favail > 0.)
		die("Insufficient atomic orbitals to occupy %lg electrons (%lg excess electrons) [per spin channel].\n", N, Favail);
	double spinFactor = (spinCount>1 && magneticMoment) ? 1. : 2.; //if unpolarized, print total occupations over both spin channels
	const char orbCode[] = "spdfgh";
	for(unsigned l=0; l<psiRadial.size(); l++) if(psiRadial[l].size())
	{	if(isRelativistic() && l>0)
		{	int l2 = 2*l;
			for(int j2=l2-1; j2<=l2+1; j2+=2)
			{	logPrintf("  %c%c (", orbCode[l], (j2<l2 ? '-' : '+'));
				for(unsigned p=0; p<psiRadial[l].size(); p++)
					if(psi2j[l][p]==j2) logPrintf(" %.2lg", atomF[l][p] * spinFactor);
				logPrintf(" )");
			}
		}
		else
		{	logPrintf("  %c (", orbCode[l]);
			for(unsigned p=0; p<psiRadial[l].size(); p++)
				logPrintf(" %.2lg", atomF[l][p] * spinFactor);
			logPrintf(" )");
		}
	}
	logPrintf("\n");
	
	//Compute atom electron density: (NOTE: assumes all orbitals are on the same grid)
	RadialFunctionR n = *psiRadial[0][0].rFunc; n.f.assign(n.r.size(), 0.); //same grid as orbitals, initialize to 0.
	for(unsigned l=0; l<psiRadial.size(); l++)
	{	for(unsigned p=0; p<psiRadial[l].size(); p++)
		{	const RadialFunctionR& psi = *(psiRadial[l][p].rFunc);
			for(unsigned i=0; i<n.r.size(); i++)
				n.f[i] += atomF[l][p] * psi.f[i] * psi.f[i];
		}
		//Augmentation for ultrasofts:
		if(Qint.size())
		{	for(unsigned p1=0; p1<VnlRadial[l].size(); p1++)
				for(unsigned p2=0; p2<=p1; p2++)
				{	QijIndex qIndex = { (int)l, (int)p1, (int)l, (int)p2, 0 };
					auto Qijl = Qradial.find(qIndex);
					if(Qijl==Qradial.end()) continue; //no augmentation for this combination
					//Collect density matrix element for this pair of projectors:
					double matrixElem = 0.;
					for(unsigned o=0; o<psiRadial[l].size(); o++)
						matrixElem += atomF[l][o]
							* dot(*psiRadial[l][o].rFunc, *VnlRadial[l][p1].rFunc)
							* dot(*psiRadial[l][o].rFunc, *VnlRadial[l][p2].rFunc);
					//Augment density:
					axpy(matrixElem * (p1==p2 ? 1. : 2.), *Qijl->second.rFunc, n);
				}
		}
	}
	//Fix normalization:
	double normFac = std::pow(e->gInfo.detR,2)/(4.*M_PI);
	for(unsigned i=0; i<n.r.size(); i++) n.f[i] *= normFac;
	//Transform density:
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	n.transform(0, dG, nGridLoc, nRadial);
}

void SpeciesInfo::getAtomPotential(RadialFunctionG& dRadial) const
{	//Get the radial density of a neutral atom:
	RadialFunctionG nRadial;
	logSuspend();
	getAtom_nRadial(0,0, nRadial, true);
	logResume();
	RadialFunctionR n = *nRadial.rFunc;
	nRadial.free();
	//Calculate cumulative charge distribution:
	std::vector<double> Qin(n.r.size());
	double Qcur = 0.;
	for(size_t i=0; i<Qin.size(); i++)
	{	double dQ = (2*n.f[i]) * (4*M_PI * n.r[i]*n.r[i] * n.dr[i]); //factor of 2 for spin
		Qcur += 0.5*dQ;
		Qin[i] = Qcur; //Trapezoidal rule
		Qcur += 0.5*dQ;
	}
	//Calculate electrostatic potential:
	RadialFunctionR d = n;
	double phiCur = Z/n.r.back();
	int iLast = int(Qin.size())-1;
	for(int i=iLast; i>=0; i--)
	{	double rInv = n.r[i] ? 1./n.r[i] : 0.;
		double dphi = Qin[i] * rInv*rInv * n.dr[i];
		if(i<iLast)
			phiCur += 0.5*dphi;
		d.f[i] = phiCur - Z*rInv;
		phiCur += 0.5*dphi;
	}
	//Pseudize potential:
	//--- calculate pseudization radius
	double Rvdw = 0.;
	for(int i=iLast; i>=0; i--)
		if(2.*n.f[i] > 0.01)
		{	Rvdw = 1.1 * n.r[i]; //DFT+D2 definition for convenience
			break;
		}
	int iVdw = std::upper_bound(n.r.begin(), n.r.end(), Rvdw) - n.r.begin();
	assert(iVdw < iLast);
	Rvdw = n.r[iVdw]; //round to grid point
	//--- pseudize electrostatic potential so that it is zero outside Rvdw
	double d0 = d.f[iVdw]; //value
	double d1 = (d.f[iVdw+1] - d.f[iVdw-1]) / log(d.r[iVdw+1]/d.r[iVdw-1]); //first derivative w.r.t x = r/Rvdw
	double d2 = (d.f[iVdw+1] + d.f[iVdw-1] - 2.*d0) / std::pow(0.5*log(d.r[iVdw+1]/d.r[iVdw-1]),2) - d1; //second derivative w.r.t x = r/Rvdw
	double a3 = 10.*d0 - 4.*d1 + 0.5*d2; //expansion x^3 (a3 + a4 x + a5 x^2) that is C2 continuous at 0 and Rvdw
	double a4 = -15.*d0 + 7.*d1 - d2;
	double a5 = 6.*d0 - 3.*d1 + 0.5*d2;
	for(int i=0; i<iVdw; i++)
	{	double x = d.r[i] / Rvdw;
		d.f[i] -= x*x*x*(a3 + x*(a4 + x*a5));
	}
	for(size_t i=iVdw; i<d.f.size(); i++)
		d.f[i] = 0.;
	//Add non-electrostatic part of Vloc:
	RadialFunctionR& Vloc = *(VlocRadial.rFunc);
	for(size_t i=0; i<d.f.size(); i++)
		d.f[i] += Vloc.f[i];
	//Transform potential:
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	d.transform(0, dG, nGridLoc, dRadial);
}
