/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <electronic/DOS.h>
#include <electronic/TetrahedralDOS.h>
#include <electronic/Everything.h>
#include <core/SphericalHarmonics.h>
#include <electronic/ColumnBundle.h>
#include <core/ScalarFieldIO.h>
#include <core/LatticeUtils.h>
#include <array>

DOS::DOS() : Etol(1e-6), Esigma(0)
{
}

void DOS::setup(const Everything& everything)
{	e = &everything;
	if(!weights.size()) //Add the default of Total Complete DOS:
		weights.push_back(Weight());
	//Check compatibility of orbital modes with pseudopotentials:
	for(const Weight& weight: weights)
		if(weight.type == Weight::Orbital)
		{	const SpeciesInfo& sp = *(e->iInfo.species[weight.specieIndex]);
			const Weight::OrbitalDesc& oDesc = weight.orbitalDesc;
			int lMax = sp.lMaxAtomicOrbitals();
			if(lMax < 0)
				die("Species %s has no atomic orbitals and cannot be used for orbital-projected density of states.\n", sp.name.c_str());
			if(oDesc.l > lMax)
				die("Angular momentum of density-of-states projection orbital %s exceeds lMax=%d of species %s.\n",
					string(oDesc).c_str(), lMax, sp.name.c_str());
			int nAOl = sp.nAtomicOrbitals(oDesc.l);
			if(oDesc.n >= unsigned(nAOl))
				die("Principal quantum number (%d) of density-of-states projection orbital %s exceeds maximum value (%d) for l=%d of species %s.\n",
					oDesc.n+1, string(oDesc).c_str(), nAOl, oDesc.l, sp.name.c_str());
		}
	Citations::add("Linear-tetrahedron sampling for density of states",
		"G. Lehmann and M. Taut, Phys. status solidi (b) 54, 469 (1972)");
}

string DOS::Weight::getDescription(const Everything& e) const
{	ostringstream oss;
	if(type==Total) oss << "Total";
	if(type==Slice || type==AtomSlice)
		oss << "(" << direction[0] << "," << direction[1] << "," << direction[2] << ") slice of half-width ";
	if(type==Sphere || type==AtomSphere) oss << "Sphere of radius ";
	if(type==Slice || type==AtomSlice || type==Sphere || type==AtomSphere)
		oss << radius << " bohr at ";
	if(type==Orbital)
		oss << orbitalDesc << " orbital at ";
	if(type==OrthoOrbital)
		oss << orbitalDesc << " orthonormalized-orbital at ";
	if(type==Sphere || type==Slice)
	{	vector3<> c;
		if(e.iInfo.coordsType == CoordsLattice)
		{	c = center;
			oss << "lattice";
		}
		else
		{	c = e.gInfo.invR * center;
			oss << "cartesian";
		}
		oss << " (" << c[0] << "," << c[1] << "," << c[2] << ")";
	}
	if(type==AtomSphere || type==AtomSlice || type==Orbital || type==OrthoOrbital)
		oss << e.iInfo.species[specieIndex]->name << " #" << (atomIndex+1);
	if(type==File)
		oss << "Weighted by '" << filename << "'";
	if(fillingMode == Occupied)
		oss << " (Occupied)";
	return oss.str();
}

//Connection between string representations of orbitals and (l,m) pairs:
const EnumStringMap<DOS::Weight::OrbitalDesc>& getOrbitalDescMap()
{	static EnumStringMap<DOS::Weight::OrbitalDesc> orbitalDescMap
	(	//Spin-less specification (allowed for all spin modes for nonrelativistic pseudopotentials)
		DOS::Weight::OrbitalDesc(0,0, 0, SpinNone), "s",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinNone), "py",
		DOS::Weight::OrbitalDesc(1, 0, 0, SpinNone), "pz",
		DOS::Weight::OrbitalDesc(1, 1, 0, SpinNone), "px",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinNone), "p",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinNone), "dxy",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinNone), "dyz",
		DOS::Weight::OrbitalDesc(2, 0, 0, SpinNone), "dz2",
		DOS::Weight::OrbitalDesc(2, 1, 0, SpinNone), "dxz",
		DOS::Weight::OrbitalDesc(2, 2, 0, SpinNone), "dx2-y2",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinNone), "d",
		DOS::Weight::OrbitalDesc(2,4, 0, SpinNone), "t2g",
		DOS::Weight::OrbitalDesc(2,5, 0, SpinNone), "eg",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinNone), "fy(3x2-y2)",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinNone), "fxyz",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinNone), "fyz2",
		DOS::Weight::OrbitalDesc(3, 0, 0, SpinNone), "fz3",
		DOS::Weight::OrbitalDesc(3, 1, 0, SpinNone), "fxz2",
		DOS::Weight::OrbitalDesc(3, 2, 0, SpinNone), "fz(x2-y2)",
		DOS::Weight::OrbitalDesc(3, 3, 0, SpinNone), "fx(x2-3y2)",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinNone), "f",
		//Up/dn specification: noncollinear modes with nonrelativistic pseudopotentials only
		DOS::Weight::OrbitalDesc(0,0, 0, SpinZ), "sUp",
		DOS::Weight::OrbitalDesc(0,0, 1, SpinZ), "sDn",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinZ), "pyUp",
		DOS::Weight::OrbitalDesc(1,-1, 1, SpinZ), "pyDn",
		DOS::Weight::OrbitalDesc(1, 0, 0, SpinZ), "pzUp",
		DOS::Weight::OrbitalDesc(1, 0, 1, SpinZ), "pzDn",
		DOS::Weight::OrbitalDesc(1, 1, 0, SpinZ), "pxUp",
		DOS::Weight::OrbitalDesc(1, 1, 1, SpinZ), "pxDn",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinZ), "pUp",
		DOS::Weight::OrbitalDesc(1,2, 1, SpinZ), "pDn",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinZ), "dxyUp",
		DOS::Weight::OrbitalDesc(2,-2, 1, SpinZ), "dxyDn",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinZ), "dyzUp",
		DOS::Weight::OrbitalDesc(2,-1, 1, SpinZ), "dyzDn",
		DOS::Weight::OrbitalDesc(2, 0, 0, SpinZ), "dz2Up",
		DOS::Weight::OrbitalDesc(2, 0, 1, SpinZ), "dz2Dn",
		DOS::Weight::OrbitalDesc(2, 1, 0, SpinZ), "dxzUp",
		DOS::Weight::OrbitalDesc(2, 1, 1, SpinZ), "dxzDn",
		DOS::Weight::OrbitalDesc(2, 2, 0, SpinZ), "dx2-y2Up",
		DOS::Weight::OrbitalDesc(2, 2, 1, SpinZ), "dx2-y2Dn",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinZ), "dUp",
		DOS::Weight::OrbitalDesc(2,3, 1, SpinZ), "dDn",
		DOS::Weight::OrbitalDesc(2,4, 0, SpinZ), "t2gUp",
		DOS::Weight::OrbitalDesc(2,4, 1, SpinZ), "t2gDn",
		DOS::Weight::OrbitalDesc(2,5, 0, SpinZ), "egUp",
		DOS::Weight::OrbitalDesc(2,5, 1, SpinZ), "egDn",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinZ), "fy(3x2-y2)Up",
		DOS::Weight::OrbitalDesc(3,-3, 1, SpinZ), "fy(3x2-y2)Dn",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinZ), "fxyzUp",
		DOS::Weight::OrbitalDesc(3,-2, 1, SpinZ), "fxyzDn",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinZ), "fyz2Up",
		DOS::Weight::OrbitalDesc(3,-1, 1, SpinZ), "fyz2Dn",
		DOS::Weight::OrbitalDesc(3, 0, 0, SpinZ), "fz3Up",
		DOS::Weight::OrbitalDesc(3, 0, 1, SpinZ), "fz3Dn",
		DOS::Weight::OrbitalDesc(3, 1, 0, SpinZ), "fxz2Up",
		DOS::Weight::OrbitalDesc(3, 1, 1, SpinZ), "fxz2Dn",
		DOS::Weight::OrbitalDesc(3, 2, 0, SpinZ), "fz(x2-y2)Up",
		DOS::Weight::OrbitalDesc(3, 2, 1, SpinZ), "fz(x2-y2)Dn",
		DOS::Weight::OrbitalDesc(3, 3, 0, SpinZ), "fx(x2-3y2)Up",
		DOS::Weight::OrbitalDesc(3, 3, 1, SpinZ), "fx(x2-3y2)Dn",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinZ), "fUp",
		DOS::Weight::OrbitalDesc(3,4, 1, SpinZ), "fDn",
		//j,mj specification (only option for spin-split channels of relativsitic pseudopotentials)
		DOS::Weight::OrbitalDesc(0,-1, 0, SpinOrbit), "s(-1/2)",
		DOS::Weight::OrbitalDesc(0,+0, 0, SpinOrbit), "s(+1/2)",
		DOS::Weight::OrbitalDesc(1,-2, 0, SpinOrbit), "p+(-3/2)",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinOrbit), "p+(-1/2)",
		DOS::Weight::OrbitalDesc(1,+0, 0, SpinOrbit), "p+(+1/2)",
		DOS::Weight::OrbitalDesc(1,+1, 0, SpinOrbit), "p+(+3/2)",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinOrbit), "p+",
		DOS::Weight::OrbitalDesc(1,+0, 1, SpinOrbit), "p-(-1/2)",
		DOS::Weight::OrbitalDesc(1,+1, 1, SpinOrbit), "p-(+1/2)",
		DOS::Weight::OrbitalDesc(1,2, 1, SpinOrbit), "p-",
		DOS::Weight::OrbitalDesc(2,-3, 0, SpinOrbit), "d+(-5/2)",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinOrbit), "d+(-3/2)",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinOrbit), "d+(-1/2)",
		DOS::Weight::OrbitalDesc(2,+0, 0, SpinOrbit), "d+(+1/2)",
		DOS::Weight::OrbitalDesc(2,+1, 0, SpinOrbit), "d+(+3/2)",
		DOS::Weight::OrbitalDesc(2,+2, 0, SpinOrbit), "d+(+5/2)",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinOrbit), "d+",
		DOS::Weight::OrbitalDesc(2,-1, 1, SpinOrbit), "d-(-3/2)",
		DOS::Weight::OrbitalDesc(2,+0, 1, SpinOrbit), "d-(-1/2)",
		DOS::Weight::OrbitalDesc(2,+1, 1, SpinOrbit), "d-(+1/2)",
		DOS::Weight::OrbitalDesc(2,+2, 1, SpinOrbit), "d-(+3/2)",
		DOS::Weight::OrbitalDesc(2,3, 1, SpinOrbit), "d-",
		DOS::Weight::OrbitalDesc(3,-4, 0, SpinOrbit), "f+(-7/2)",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinOrbit), "f+(-5/2)",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinOrbit), "f+(-3/2)",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinOrbit), "f+(-1/2)",
		DOS::Weight::OrbitalDesc(3,+0, 0, SpinOrbit), "f+(+1/2)",
		DOS::Weight::OrbitalDesc(3,+1, 0, SpinOrbit), "f+(+3/2)",
		DOS::Weight::OrbitalDesc(3,+2, 0, SpinOrbit), "f+(+5/2)",
		DOS::Weight::OrbitalDesc(3,+3, 0, SpinOrbit), "f+(+7/2)",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinOrbit), "f+",
		DOS::Weight::OrbitalDesc(3,-2, 1, SpinOrbit), "f-(-5/2)",
		DOS::Weight::OrbitalDesc(3,-1, 1, SpinOrbit), "f-(-3/2)",
		DOS::Weight::OrbitalDesc(3,+0, 1, SpinOrbit), "f-(-1/2)",
		DOS::Weight::OrbitalDesc(3,+1, 1, SpinOrbit), "f-(+1/2)",
		DOS::Weight::OrbitalDesc(3,+2, 1, SpinOrbit), "f-(+3/2)",
		DOS::Weight::OrbitalDesc(3,+3, 1, SpinOrbit), "f-(+5/2)",
		DOS::Weight::OrbitalDesc(3,4, 1, SpinOrbit), "f-"
	);
	return orbitalDescMap;
}

void DOS::Weight::OrbitalDesc::parse(string desc)
{	//SPlit into principal quantum number and orbital code:
	size_t orbStart = desc.find_first_not_of("0123456789");
	if(orbStart == string::npos)
		throw "Orbital description '" + desc + "' does not have an orbital code.";
	string orbCode = desc.substr(orbStart);
	//Convert orbital code:
	if(!getOrbitalDescMap().getEnum(orbCode.c_str(), *this))
		throw "'" + orbCode + "' is not a valid orbital code for an orbital description.";
	//Extract the principal quantum number:
	if(orbStart)
	{	int nPlus1 = atoi(desc.substr(0, orbStart).c_str());
		if(nPlus1<=0) throw string("Principal quantum number in orbital description must be a positive integer");
		n = nPlus1 - 1;
	}
}

DOS::Weight::OrbitalDesc::operator string() const
{	ostringstream oss;
	if(n>0) oss << (n+1); //optional pseudo-principal quantum number 
	oss << getOrbitalDescMap().getString(*this);
	return oss.str();
}

bool DOS::Weight::OrbitalDesc::operator<(const DOS::Weight::OrbitalDesc& other) const
{	if(l != other.l) return l < other.l;
	if(m != other.m) return m < other.m;
	if(s != other.s) return s < other.s;
	if(n != other.n) return n < other.n;
	return spinType < other.spinType;
}

//Thread function for setting fourier transform of slice of half-width R centered at r0 parallel to lattice-plane d:
inline void sliceWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
	complex* w, const vector3<>& r0, double R, const vector3<int>& d)
{
	double planeSpacing = (2*M_PI) / sqrt(GGT.metric_length_squared(gcdReduce(d))); //lattice-plane spacing
	double prefac = (2*R) / planeSpacing;
	int dSq = d.length_squared();
	THREAD_halfGspaceLoop
	(	int iGsq = iG.length_squared();
		int iGdotd = dot(iG, d);
		if(iGsq * dSq == iGdotd * iGdotd) // => iG || d (Cauchy-Schwarz)
		{	double GR = R * sqrt(GGT.metric_length_squared(iG));
			w[i] = prefac * cis(-(2.*M_PI)*dot(iG,r0)) * bessel_jl(0,GR);
		}
		else w[i] = 0.;
	)
}

//Thread function for setting fourier transform of sphere of radius R centered at r0:
inline void sphereWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
	complex* w, const vector3<>& r0, double R, double cellVolume)
{
	double prefac = (4./3)*M_PI * pow(R,3) / cellVolume;
	THREAD_halfGspaceLoop
	(	double GR = R * sqrt(GGT.metric_length_squared(iG));
		w[i] = prefac * cis(-(2.*M_PI)*dot(iG,r0)) * (bessel_jl(0,GR) + bessel_jl(2,GR));
	)
}


void DOS::dump()
{
	const ElecInfo& eInfo = e->eInfo;
	int nSpins = eInfo.nSpins();
	const Supercell& supercell = *(e->coulombParams.supercell);
	
	//Initialize evaluator:
	std::vector<int> iReduced(supercell.kmesh.size());
	for(size_t ik=0; ik<iReduced.size(); ik++)
		iReduced[ik] = supercell.kmeshTransform[ik].iReduced;
	TetrahedralDOS eval(supercell.kmesh, iReduced, e->gInfo.R, supercell.super,
		nSpins, eInfo.nBands, weights.size(), eInfo.qWeightSum/nSpins);
	
	//Compute the weight functions on the real-space grid for weighted-density modes:
	std::vector<ScalarFieldArray> weightFuncs(weights.size());
	bool needDensity = false; //check if any of the modes need the density
	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
	{	Weight& weight = weights[iWeight];
		ScalarFieldArray& weightFunc = weightFuncs[iWeight];
		switch(weight.type)
		{	case Weight::Slice:
			case Weight::Sphere:
			case Weight::AtomSlice:
			case Weight::AtomSphere:
			{	//Get the center from atom position for Atom* modes:
				if(weight.type == Weight::AtomSlice || weight.type == Weight::AtomSphere)
					weight.center = e->iInfo.species[weight.specieIndex]->atpos[weight.atomIndex];
				//Compute weight function in Fourier space:
				const GridInfo& gInfo = e->gInfo;
				ScalarFieldTilde weightFuncTilde(ScalarFieldTildeData::alloc(gInfo));
				complex* wData = weightFuncTilde->data();
				if(weight.type == Weight::Slice || weight.type == Weight::AtomSlice)
					threadLaunch(sliceWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, weight.direction);
				else // (weight.type == Weight::Sphere || weight.type == Weight::Atomphere)
					threadLaunch(sphereWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, fabs(gInfo.detR));
				//Store weight function in real space:
				weightFunc.resize(1);
				weightFunc[0] = I(weightFuncTilde);
				needDensity = true;
				break;
			}
			case Weight::File:
			{	nullToZero(weightFunc, e->gInfo, 1);
				loadRawBinary(weightFunc[0], weight.filename.c_str());
				needDensity = true;
				break;
			}
			default: //Total and orbital modes do not weight the density:
			{	break;
			}
		}
		if(eInfo.spinType == SpinZ && weightFunc.size())
		{	weightFunc.resize(eInfo.nDensities);
			weightFunc[1] = weightFunc[0];
		}
		if(eInfo.spinType == SpinVector)
		{	if(weight.Mhat.length() && weight.type==Weight::Total)
			{	//treat spin-projected total DOS as a density-weighted DOS with spatially uniform weight
				nullToZero(weightFunc, e->gInfo, 1);
				weightFunc[0] += 1.;
				needDensity = true;
			}
			if(weightFunc.size()) //convert to the appropriate 4-array depending on magnetization
			{	ScalarField w; std::swap(w, weightFunc[0]);
				weightFunc.resize(eInfo.nDensities);
				const vector3<>& M = weight.Mhat;
				weightFunc[0] = w * (1. + M[2]); //UpUp
				weightFunc[1] = w * (1. - M[2]); //DnDn
				weightFunc[2] = w * (+2.*M[0]); //Re(UpDn)
				weightFunc[3] = w * (-2.*M[1]); //Im(UpDn)
				if(M.length_squared()) weightFunc *= 0.5; //count only one spin
			}
		}
	}
	//Compute the overlap of the density for each states and band, with each of the density weight functions:
	if(needDensity)
	{	//Direct contributions from bands:
		for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		{	for(int iBand=0; iBand<eInfo.nBands; iBand++)
			{	//Compute the density for this state and band:
				diagMatrix F(1, 1.); //compute density with filling=1; incorporate fillings later per weight function if required
				ColumnBundle C = e->eVars.C[iState].getSub(iBand, iBand+1);
				ScalarFieldArray n = diagouterI(F, C, eInfo.nDensities, &e->gInfo);
				//Compute the weights:
				for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
					if(weightFuncs[iWeight].size())
						eval.w(iWeight, iState, iBand) = e->gInfo.dV * dot(weightFuncs[iWeight], n);
			}
		}
		//Ultrasoft augmentation:
		for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
			if(weightFuncs[iWeight].size())
			{	//Project weight functions from grid to atom-centered spherical functions:
				e->iInfo.augmentDensityGridGrad(weightFuncs[iWeight]);
				//Propagate from spherical functions to bands:
				for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
				{	const std::vector<matrix>& VdagCq = e->eVars.VdagC[iState];
					std::vector<matrix> wVdagCq(e->iInfo.species.size());
					const QuantumNumber& qnum = e->eInfo.qnums[iState];
					e->iInfo.augmentDensitySphericalGrad(qnum, VdagCq, wVdagCq);
					//Sum diagonal contributions per band over species
					diagMatrix wAug(e->eInfo.nBands, 0.);
					for(size_t sp=0; sp<VdagCq.size(); sp++)
						if(wVdagCq[sp]) wAug += diag(dagger(VdagCq[sp]) * wVdagCq[sp]);
					//Augment with appropriate prefactors:
					for(int iBand=0; iBand<eInfo.nBands; iBand++)
						eval.w(iWeight, iState, iBand) += e->gInfo.dV * wAug[iBand];
				}
			}
	}
	
	//Check which orbital modes are needed:
	bool needOrbitals = false;
	bool needOrthoOrbitals = false;
	for(const Weight& weight: weights)
	{	if(weight.type == Weight::Orbital) needOrbitals = true;
		if(weight.type == Weight::OrthoOrbital) needOrthoOrbitals = true;
	}
	
	//Compute projections for orbital mode:
	if(needOrbitals || needOrthoOrbitals)
	{	for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		{	const ColumnBundle& C = e->eVars.C[iState];
			//Count atomic orbitals:
			int nOrbitals = 0;
			std::vector<int> spOffset(e->iInfo.species.size()); //species offset into atomic orbitals list
			for(unsigned sp=0; sp<e->iInfo.species.size(); sp++)
			{	spOffset[sp] = nOrbitals;
				nOrbitals += e->iInfo.species[sp]->nAtomicOrbitals();
			}
			//Atomic-orbital projections if needed:
			matrix CdagOpsi;
			if(needOrbitals)
				CdagOpsi = C ^ e->iInfo.getAtomicOrbitals(iState, true); //directly calculate overlap using Opsi
			//Ortho-orbital projections if needed:
			matrix CdagOpsiOrtho;
			if(needOrthoOrbitals)
			{	ColumnBundle psi = e->iInfo.getAtomicOrbitals(iState, false);
				//Orthogonalize:
				ColumnBundle Opsi = O(psi);
				matrix orthoMat = invsqrt(psi ^ Opsi); //orthonormalizing matrix
				CdagOpsiOrtho = (C ^ Opsi) * orthoMat;
			}
			for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
			{	const Weight& weight = weights[iWeight];
				if(weight.type == Weight::Orbital || weight.type == Weight::OrthoOrbital)
				{	const Weight::OrbitalDesc& oDesc = weight.orbitalDesc;
					const matrix& CdagOpsiRel = weight.type==Weight::Orbital ? CdagOpsi : CdagOpsiOrtho;
					int l = oDesc.l;
					//Check relativistic psp compatibility:
					if(e->iInfo.species[weight.specieIndex]->isRelativistic() && oDesc.l)
					{	if(!(oDesc.spinType==SpinOrbit || (oDesc.spinType==SpinNone && oDesc.m==l+1)))
							die("Individual orbital projections for l>0 in relativsitic pseudopotentials must use the j,mj specification (eg. 'p+(+1/2)').\n");
					}
					else
					{	if(oDesc.spinType==SpinOrbit)
							die("Only l>0 projections in relativistic pseudopotentials may use the j,mj specification.\n");
					}
					//Initialize to zero (eval has 1 by default):
					for(int iBand=0; iBand<eInfo.nBands; iBand++)
						eval.w(iWeight, iState, iBand) = 0.;
					//Accumulate weights:
					int sStart, sStop;
					if(oDesc.spinType==SpinNone)
					{	sStart=0;
						sStop=eInfo.spinorLength();
					}
					else
					{	sStart=oDesc.s;
						sStop=oDesc.s+1;
					}
					for(int s=sStart; s<sStop; s++)
					{	std::vector<int> mArr; //set of m to include
						if(oDesc.m==l+1) //all orbitals
						{	int mMin = -l;
							int mMax = +l;
							if(e->iInfo.species[weight.specieIndex]->isRelativistic())
								mMin -= (s ? -1 : +1);
							for(int m=mMin; m<=mMax; m++)
								mArr.push_back(m);
						}
						else if(l==2 && oDesc.m==l+2) //t2g set of d orbitals
						{	mArr.push_back(-2); //dxy
							mArr.push_back(-1); //dyz
							mArr.push_back(+1); //dxz
						}
						else if(l==2 && oDesc.m==l+3) //eg set of d orbitals
						{	mArr.push_back(0); //dz2
							mArr.push_back(2); //dx2-y2
						}
						else //single orbital mode
						{	mArr.push_back(oDesc.m);
						}
						for(int m: mArr)
						{	int iCol = spOffset[weight.specieIndex] + e->iInfo.species[weight.specieIndex]->atomicOrbitalOffset(weight.atomIndex, oDesc.n, l, m, s);
							for(int iBand=0; iBand<eInfo.nBands; iBand++)
								eval.w(iWeight, iState, iBand) += (CdagOpsiRel.data()[CdagOpsiRel.index(iBand, iCol)]).norm();
						}
					}
				}
			}
		}
	}
	
	//Read override eigenvalues if any:
	const std::vector<diagMatrix>* eigs = &e->eVars.Hsub_eigs;
	std::vector<diagMatrix> eigsOverride;
	if(eigsFilename.length())
	{	//Read over-ride eigenvalues from specified file:
		logPrintf("Reading over-ride eigenvalues from file %s ... ", eigsFilename.c_str()); logFlush();
		eigsOverride.resize(eInfo.nStates);
		eInfo.read(eigsOverride, eigsFilename.c_str());
		eigs = &eigsOverride;
		logPrintf("done.\n");
	}
	
	//Apply the filling weights and set the eigenvalues:
	for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		for(int iBand=0; iBand<eInfo.nBands; iBand++)
		{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
				if(weights[iWeight].fillingMode == Weight::Occupied)
					eval.w(iWeight, iState, iBand) *= e->eVars.F[iState][iBand];
			eval.e(iState, iBand) = (*eigs)[iState][iBand];
		}
		
	//Synchronize eigenvalues and weights between processes:
	if(mpiWorld->nProcesses()>1)
	{	if(mpiWorld->isHead())
		{	for(int iSrc=1; iSrc<mpiWorld->nProcesses(); iSrc++)
			{	int qStart = eInfo.qStartOther(iSrc);
				int qStop = eInfo.qStopOther(iSrc);
				std::vector<double> message((qStop-qStart)*eInfo.nBands*(weights.size()+1));
				mpiWorld->recv(message.data(), message.size(), iSrc, 0);
				const double* messagePtr = message.data();
				for(int iState=qStart; iState<qStop; iState++)
					for(int iBand=0; iBand<eInfo.nBands; iBand++)
					{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
							eval.w(iWeight, iState, iBand) = *(messagePtr++);
						eval.e(iState, iBand) = *(messagePtr++);
					}
			}
		}
		else
		{	//Pack data into a single message:
			std::vector<double> message((eInfo.qStop-eInfo.qStart)*eInfo.nBands*(weights.size()+1));
			double* messagePtr = message.data();
			for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
				for(int iBand=0; iBand<eInfo.nBands; iBand++)
				{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
						*(messagePtr++) = eval.w(iWeight, iState, iBand);
					*(messagePtr++) = eval.e(iState, iBand);
				}
			mpiWorld->send(message.data(), message.size(), 0, 0);
		}
	}
	if(!mpiWorld->isHead()) return;
	
	//Compute and print density of states (head only):
	string header = "\"Energy\"";
	for(const Weight& weight: weights)
		header += ("\t\"" + weight.getDescription(*e) + "\"");
	eval.weldEigenvalues(Etol);
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{	TetrahedralDOS::Lspline dos = eval.getDOS(iSpin, Etol);
		if(Esigma>0.) dos = eval.gaussSmooth(dos, Esigma); //apply Gauss smoothing if requested
		eval.printDOS(dos, e->dump.getFilename(nSpins==1 ? "dos" : (iSpin==0 ? "dosUp" : "dosDn")), header);
	}
}
