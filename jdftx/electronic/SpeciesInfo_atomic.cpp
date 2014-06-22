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
#include <core/LatticeUtils.h>

//------- SpeciesInfo functions related to atomic orbitals -------


//RadialFunctionR operators implemented in SpeciesInfo_readUSPP
double dot(const RadialFunctionR& X, const RadialFunctionR& Y);
void axpy(double alpha, const RadialFunctionR& X, RadialFunctionR& Y);


void SpeciesInfo::accumulateAtomicDensity(DataGptrCollection& nTilde) const
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
	//Scratch space for atom positions on GPU:
	#ifdef GPU_ENABLED
	vector3<>* atposCurPref;
	cudaMalloc(&atposCurPref, sizeof(vector3<>)*atpos.size());
	#endif
	//For each magnetization magnitude:
	for(const MomentMagnitude& Mmag: Mmags)
	{	if(e->eInfo.nDensities == 1) assert(!Mmag.Mlength);
		//Get the radial spin-densities in diagonal basis:
		std::vector<RadialFunctionG> nRadial(Mmag.Mlength ? 2 : 1);
		for(unsigned s=0; s<nRadial.size(); s++)
			getAtom_nRadial(s, Mmag.Mlength, nRadial[s]);
		//Collect contributions for each direction with this magnitude:
		for(const MomentDirection& Mdir: Mmag.Mdirs)
		{	//Compute structure factor for atoms with current magnetization vector:
			#ifdef GPU_ENABLED
			cudaMemcpy(atposCurPref, Mdir.atpos.data(), Mdir.atpos.size()*sizeof(vector3<>), cudaMemcpyHostToDevice);
			#else
			const vector3<>* atposCurPref = Mdir.atpos.data();
			#endif
			DataGptr SG; nullToZero(SG, e->gInfo);
			callPref(getSG)(e->gInfo.S, Mdir.atpos.size(), atposCurPref, 1./e->gInfo.detR, SG->dataPref());
			//Spin-densities in diagonal basis
			std::vector<DataGptr> nDiag;
			for(const RadialFunctionG& nRad: nRadial)
				nDiag.push_back(nRad * SG);
			//Accumulate contributions:
			if(nDiag.size()==1)
			{	if(nTilde.size()==1)
					nTilde[0] += nDiag[0];
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
	//Cleanup:
	#ifdef GPU_ENABLED
	cudaFree(atposCurPref);
	#endif
}

//Set atomic orbitals in column bundle from radial functions (almost same operation as setting Vnl)
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& Y, int colOffset) const
{	if(!atpos.size()) return;
	assert(Y.basis); assert(Y.qnum);
	//Check sizes:
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	if(nSpinCopies>1) assert(Y.isSpinor()); //can have multiple spinor copies only in spinor mode
	int colMax = colOffset;
	for(int l=0; l<int(psiRadial.size()); l++)
		colMax += psiRadial[l].size() * (2*l+1) * atpos.size() * nSpinCopies;
	assert(colMax <= Y.nCols());
	//Set orbitals and associated info if requested:
	const Basis& basis = *Y.basis;
	int iCol = colOffset / nSpinCopies; //current column (counting equal and opposite spinors as one)
	for(int l=0; l<int(psiRadial.size()); l++)
		for(unsigned p=0; p<psiRadial[l].size(); p++)
		{	for(int m=-l; m<=l; m++)
			{	//Set atomic orbitals for all atoms at specified (n,l,m):
				size_t atomStride = Y.colLength() * nSpinCopies;
				size_t offs = iCol * atomStride;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G, atposPref, psiRadial[l][p], Y.dataPref()+offs);
				if(nSpinCopies>1) //make copy for other spin and zero the minor component
				{	complex* dataPtr = Y.dataPref()+offs;
					for(size_t a=0; a<atpos.size(); a++)
					{	callPref(eblas_zero)(2*basis.nbasis, dataPtr+basis.nbasis);
						callPref(eblas_copy)(dataPtr+3*basis.nbasis, dataPtr, basis.nbasis);
						dataPtr += atomStride;
					}
				}
				iCol += atpos.size();
			}
		}
}
int SpeciesInfo::nAtomicOrbitals() const
{	int nOrbitals = 0;
	for(int l=0; l<int(psiRadial.size()); l++)
		nOrbitals += (2*l+1)*psiRadial[l].size();
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	return nOrbitals * atpos.size() * nSpinCopies;
}
int SpeciesInfo::lMaxAtomicOrbitals() const
{	return int(psiRadial.size()) - 1;
}
int SpeciesInfo::nAtomicOrbitals(int l) const
{	assert(l >= 0);
	if(unsigned(l) >= psiRadial.size()) return -1; //signals end of l
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	return psiRadial[l].size() * nSpinCopies;
}
int SpeciesInfo::atomicOrbitalOffset(unsigned int iAtom, unsigned int n, int l, int m) const
{	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(n < psiRadial[l].size());
	assert(m >= -l); assert(m <= l);
	int iProj = l + m; //#projectors before this one at current l,n
	for(int L=0; L<=l; L++) //#projectors from previous l,n:
		iProj += (L==l ? n : psiRadial[l].size()) * (2*L+1);
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	return (iProj * atpos.size() + iAtom) * nSpinCopies;
}

void SpeciesInfo::setOpsi(ColumnBundle& Opsi, unsigned n, int l) const
{	if(!atpos.size()) return;
	int nSpinCopies = 2/e->eInfo.qWeightSum;
	assert(Opsi.basis); assert(Opsi.qnum);
	assert((2*l+1)*int(atpos.size())*nSpinCopies <= Opsi.nCols());
	if(nSpinCopies>1) assert(Opsi.isSpinor()); //can have multiple spinor copies only in spinor mode
	const Basis& basis = *Opsi.basis;
	int iCol = 0; //current column
	for(int m=-l; m<=l; m++)
	{	//Set atomic orbitals for all atoms at specified (n,l,m):
		size_t atomStride = Opsi.colLength() * nSpinCopies;
		size_t offs = iCol * atomStride;
		callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Opsi.qnum->k, basis.iGarrPref, e->gInfo.G, atposPref, OpsiRadial->at(l)[n], Opsi.dataPref()+offs);
		if(nSpinCopies>1) //make copy for other spin
		{	complex* dataPtr = Opsi.dataPref()+offs;
			for(size_t a=0; a<atpos.size(); a++)
			{	callPref(eblas_zero)(2*basis.nbasis, dataPtr+basis.nbasis);
				callPref(eblas_copy)(dataPtr+3*basis.nbasis, dataPtr, basis.nbasis);
				dataPtr += atomStride;
			}
		}
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
	int spinCount = (e->eInfo.nDensities==1 ? 1 : 2);
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
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	n.transform(0, dG, nGridLoc, nRadial);
}
