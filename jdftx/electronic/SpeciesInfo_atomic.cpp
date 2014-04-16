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
#include <electronic/operators.h>

//------- SpeciesInfo functions related to atomic orbitals -------


//RadialFunctionR operators implemented in SpeciesInfo_readUSPP
double dot(const RadialFunctionR& X, const RadialFunctionR& Y);
void axpy(double alpha, const RadialFunctionR& X, RadialFunctionR& Y);


void SpeciesInfo::accumulateAtomicDensity(DataGptrCollection& nTilde) const
{	//Collect list of distinct magnetizations and corresponding atom positions:
	std::map<double, std::vector< vector3<> > > Matpos;
	if(initialMagneticMoments.size())
	{	for(unsigned a=0; a<atpos.size(); a++)
			Matpos[initialMagneticMoments[a]].push_back(atpos[a]);
	}
	else Matpos[0.] = atpos;
	//Scratch space for atom positions on GPU:
	#ifdef GPU_ENABLED
	vector3<>* atposCurPref;
	cudaMalloc(&atposCurPref, sizeof(vector3<>)*atpos.size());
	#endif
	//For each magnetization:
	for(auto mapEntry: Matpos)
	{	double M = mapEntry.first; if(e->eInfo.spinType == SpinNone) assert(!M);
		const std::vector< vector3<> >& atposCur = mapEntry.second;
		#ifdef GPU_ENABLED
		cudaMemcpy(atposCurPref, atposCur.data(), atposCur.size()*sizeof(vector3<>), cudaMemcpyHostToDevice);
		#else
		const vector3<>* atposCurPref = atposCur.data();
		#endif
		//Compute structure factor for atoms with current magnetization:
		DataGptr SG; nullToZero(SG, e->gInfo);
		callPref(getSG)(e->gInfo.S, atposCur.size(), atposCurPref, 1./e->gInfo.detR, SG->dataPref());
		//Collect contributions:
		for(int s=0; s<(M ? 2 : 1); s++)
		{	RadialFunctionG nRadial;
			getAtom_nRadial(s, M, nRadial);
			DataGptr nTilde_s = nRadial * SG; //contribution per spin 
			nRadial.free();
			nTilde[s] += nTilde_s;
			if(!M) nTilde.back() += nTilde_s; //same contribution to both spins (irrespective of spinType)
		}
	}
	//Cleanup:
	#ifdef GPU_ENABLED
	cudaFree(atposCurPref);
	#endif
}

//Set atomic orbitals in column bundle from radial functions (almost same operation as setting Vnl)
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& Y, int colOffset, std::vector<AtomConfig>* atomConfig) const
{	if(!atpos.size()) return;
	assert(Y.basis); assert(Y.qnum);
	//Check sizes:
	int colMax = colOffset; for(int l=0; l<int(psiRadial.size()); l++) colMax += psiRadial[l].size() * (2*l+1) * atpos.size();
	assert(colMax <= Y.nCols());
	if(atomConfig && atomConfig->size()<size_t(colMax)) atomConfig->resize(colMax);
	//Set orbitals and associated info if requested:
	const Basis& basis = *Y.basis;
	int iCol = colOffset; //current column
	for(int l=0; l<int(psiRadial.size()); l++)
		for(unsigned p=0; p<psiRadial[l].size(); p++)
		{	for(int m=-l; m<=l; m++)
			{	//Set atomic orbitals for all atoms at specified (n,l,m):
				size_t offs = iCol * basis.nbasis;
				size_t atomStride = basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G, atposPref, psiRadial[l][p], Y.dataPref()+offs);
				if(atomConfig)
					for(unsigned a=0; a<atpos.size(); a++)
					{	AtomConfig& ac = atomConfig->at(iCol + a);
						ac.n = p;
						ac.l = l;
						ac.iAtom = a;
						//ac.F = atomF[l][p];
					}
				iCol += atpos.size();
			}
		}
}
void SpeciesInfo::setAtomicOrbital(ColumnBundle& Y, int col, unsigned iAtom, unsigned n, int l, int m) const
{	//Check inputs:
	assert(Y.basis); assert(Y.qnum);
	assert(col >= 0); assert(col < Y.nCols());
	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(n < psiRadial[l].size());
	assert(m >= -l); assert(m <= l);
	//Call Vnl() to set the column:
	const Basis& basis = *Y.basis;
	size_t offs = col * basis.nbasis;
	size_t atomStride = basis.nbasis;
	callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G, atposPref+iAtom, psiRadial[l][n], Y.dataPref()+offs);
}
int SpeciesInfo::nAtomicOrbitals() const
{	int nOrbitals = 0;
	for(int l=0; l<int(psiRadial.size()); l++)
		nOrbitals += (2*l+1)*psiRadial[l].size();
	return nOrbitals * atpos.size();
}
int SpeciesInfo::lMaxAtomicOrbitals() const
{	return int(psiRadial.size()) - 1;
}
int SpeciesInfo::nAtomicOrbitals(int l) const
{	assert(l >= 0);
	if(unsigned(l) >= psiRadial.size()) return -1; //signals end of l
	return psiRadial[l].size();
}
int SpeciesInfo::atomicOrbitalOffset(unsigned int iAtom, unsigned int n, int l, int m) const
{	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(n < psiRadial[l].size());
	assert(m >= -l); assert(m <= l);
	int iProj = l + m; //#projectors before this one at current l,n
	for(int L=0; L<=l; L++) //#projectors from previous l,n:
		iProj += (L==l ? n : psiRadial[l].size()) * (2*L+1);
	return iProj * atpos.size() + iAtom;
}

void SpeciesInfo::setOpsi(ColumnBundle& Opsi, unsigned n, int l) const
{	if(!atpos.size()) return;
	assert(Opsi.basis); assert(Opsi.qnum);
	assert((2*l+1)*int(atpos.size()) <= Opsi.nCols());
	const Basis& basis = *Opsi.basis;
	int iCol = 0; //current column
	for(int m=-l; m<=l; m++)
	{	//Set atomic orbitals for all atoms at specified (n,l,m):
		size_t offs = iCol * basis.nbasis;
		size_t atomStride = basis.nbasis;
		callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Opsi.qnum->k, basis.iGarrPref, e->gInfo.G, atposPref, OpsiRadial->at(l)[n], Opsi.dataPref()+offs);
		iCol += atpos.size();
	}
}


void SpeciesInfo::estimateAtomEigs()
{	if(!psiRadial.size()) return; //no orbitals

	//Compute eigenvalues if unavailable
	if(!atomEigs.size())
	{	std::map<int,int> invalidPsis; //some FHI pseudodpotentials store unbound projectors which need to be discarded for LCAO
		logPrintf("  Approximate pseudo-atom eigenvalues: ");
		atomEigs.resize(psiRadial.size());
		for(unsigned l=0; l<psiRadial.size(); l++)
		{	logPrintf(" l=%d: (", l);
			for(unsigned p=0; p<psiRadial[l].size(); p++)
			{	const RadialFunctionR& psi = *(psiRadial[l][p].rFunc);
				//Find points deep in the tail of the wavefunction, but still far above roundoff limit:
				double r[2], e[2];
				double psiThresh[2] = { 3e-7, 3e-6 };
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
					logPrintf(" %.2lg", eCorrected);
				}
				else
				{	psiRadial[l][p].free();
					psiRadial[l].erase(psiRadial[l].begin()+p);
					p--;
					invalidPsis[l]++;
				}
			}
			logPrintf(" )");
		}
		logPrintf("\n");
		//Report removal of invalid psi's:
		for(auto invalidPsi: invalidPsis)
			logPrintf("  WARNING: Discarded %d l=%d unbound projectors from atomic orbital set.\n", invalidPsi.second, invalidPsi.first);
	}
}

void SpeciesInfo::getAtom_nRadial(int spin, double magneticMoment, RadialFunctionG& nRadial) const
{
	int spinCount = (e->eInfo.spinType==SpinNone ? 1 : 2);
	assert(spin >= 0); assert(spin < spinCount);
	
	//Determine occupations:
	std::map< double, std::pair<unsigned,unsigned> > eigMap; //map from eigenvalues to (l,p)
	std::vector<std::vector<double> > atomF(psiRadial.size());
	for(unsigned l=0; l<psiRadial.size(); l++)
	{	atomF[l].resize(psiRadial[l].size());
		for(unsigned p=0; p<psiRadial[l].size(); p++)
			eigMap[atomEigs[l][p]] = std::make_pair(l,p);
	}
	if(spinCount>1 && magneticMoment)
		logPrintf("%s (M=%lg) pseudo-atom spin=%+d occupations: ", name.c_str(), magneticMoment, 1-2*spin);
	else
		logPrintf("%s pseudo-atom occupations: ", name.c_str());
	double Nvalence = Z - initialOxidationState;
	double N = 0.5*(Nvalence + (1-2*spin)*magneticMoment); //total electrons to fill
	if(N < 0.)
		die("Magnetic moment (%lg) exceeds pseudo-atom valence electron count (%lg) [per spin channel].\n", magneticMoment, Nvalence);
	double Favail = N; //electrons yet to be filled
	for(auto eigEntry: eigMap) //in ascending order of eigenvalues
	{	unsigned l = eigEntry.second.first;
		unsigned p = eigEntry.second.second;
		double capacity = (2*l+1);
		atomF[l][p] = std::min(Favail, capacity);
		Favail -= atomF[l][p];
	}
	if(Favail > 0.)
		die("Insufficient atomic orbitals to occupy %lg electrons (%lg excess electrons) [per spin channel].\n", N, Favail);
	double spinFactor = (spinCount>1 && magneticMoment) ? 1. : 2.; //if unpolarized, print total occupations over both spin channels
	for(unsigned l=0; l<psiRadial.size(); l++) if(psiRadial[l].size())
	{	logPrintf(" l=%d: (", l);
		for(unsigned p=0; p<psiRadial[l].size(); p++)
			logPrintf(" %.2lg", atomF[l][p] * spinFactor);
		logPrintf(" )");
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
	double normFac = 1./(4.*M_PI);
	for(unsigned i=0; i<n.r.size(); i++) n.f[i] *= normFac;
	//Transform density:
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dGloc))+5;
	n.transform(0, dGloc, nGridLoc, nRadial);
}
