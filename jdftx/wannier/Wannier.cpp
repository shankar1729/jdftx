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

#include <wannier/Wannier.h>
#include <wannier/WannierMinimizerFD.h>
#include <wannier/WannierMinimizerRS.h>

Wannier::Wannier() : needAtomicOrbitals(false), addAtomicOrbitals(false), ignoreSemiCore(true),
	localizationMeasure(LM_FiniteDifference),
	bStart(0), eOuterMin(-DBL_MAX), eOuterMax(+DBL_MAX), eInnerMin(-DBL_MAX), eInnerMax(+DBL_MAX),
	outerWindow(false), innerWindow(false), nFrozen(0),
	saveWfns(false), saveWfnsRealSpace(false), saveMomenta(false), saveSpin(false), saveRP(false),
	zFieldMag(0.),
	z0(0.), zH(0.), zSigma(0.),
	loadRotations(false), numericalOrbitalsOffset(0.5,0.5,0.5), rSmooth(1.),
	spinMode(SpinAll), polar(false)
{
}

void Wannier::setup(const Everything& everything)
{	e = &everything;
	logPrintf("\n---------- Initializing Wannier Function solver ----------\n");

	//Initialize minimization parameters:
	minParams.fpLog = globalLog;
	minParams.linePrefix = "WannierMinimize: ";
	minParams.energyLabel = "Omega";

	//Add atomic orbitals automatically, if requested:
	if(addAtomicOrbitals)
	{	logPrintf("\nAdding atomic orbitals as trial orbitals (%s semicore orbitals):\n",
			ignoreSemiCore ? "ignoring" : "including");
		for(int iSp=0; iSp<int(e->iInfo.species.size()); iSp++)
		{	const SpeciesInfo& sp = *(e->iInfo.species[iSp]);
			//Prepare template of atomic orbitals for one atom:
			DOS::Weight::OrbitalDesc od;
			od.spinType = (e->eInfo.spinType == SpinNone)
				? SpinNone
				: (sp.isRelativistic() ? SpinOrbit : SpinZ);
			std::vector<DOS::Weight::OrbitalDesc> orbitalDescs;
			for(od.l=0; od.l<=sp.lMaxAtomicOrbitals(); od.l++)
			{	int nMax = sp.nAtomicOrbitals(od.l) - 1; //max pseudo-principal quantum number
				int nMin = (ignoreSemiCore ? std::max(nMax, 0) : 0);
				for(od.n=nMin; int(od.n)<=nMax; od.n++)
				{	switch(od.spinType)
					{	case SpinNone:
						case SpinZ:
						{	int nSpins = (od.spinType == SpinNone) ? 1 : 2;
							for(od.m=-od.l; od.m<=od.l; od.m++)
								for(od.s=0; od.s<nSpins; od.s++)
									orbitalDescs.push_back(od);
							break;
						}
						case SpinOrbit:
						{	for(od.s=1; od.s>=0; od.s--) //j=l-1/2 and j=l+1/2
							{	int mStart = (od.s ? 1 : -1) - od.l;
								for(od.m=mStart; od.m<=od.l; od.m++) //mj - 1/2
									orbitalDescs.push_back(od);
							}
							break;
						}
						case SpinVector:; //Unused, but here to avoid compiler warning
					}
				}
			}
			//Add above template for each atom of species:
			TrialOrbital trialOrbital;
			trialOrbital.push_back(AtomicOrbital()); //single orbital projection
			AtomicOrbital& ao = trialOrbital[0];
			ao.sp = iSp;
			for(ao.atom=0; ao.atom<int(sp.atpos.size()); ao.atom++)
				for(DOS::Weight::OrbitalDesc od: orbitalDescs)
				{	ao.orbitalDesc = od;
					trialOrbitals.push_back(trialOrbital);
				}
			logPrintf("  Added %lu orbitals each for %lu %s atoms.\n",
				orbitalDescs.size(), sp.atpos.size(), sp.name.c_str());
		}
		logPrintf("\n");
		needAtomicOrbitals = true;
	}

	//Check window settings:
	nCenters = trialOrbitals.size() + nFrozen;
	if(nFrozen)
	{	if(outerWindow && !innerWindow)
			die("Frozen centers require an inner window, if an outer window is specified.\n");
	}
	if(innerWindow)
	{	if(!outerWindow) die("Inner window requires that an outer window be specified.\n");
		if(eInnerMin<eOuterMin || eInnerMax>eOuterMax)
			die("Inner window must lie entirely within the outer window.\n");
	}
	if(!outerWindow) //fixed bands
	{	int bStop = bStart + nCenters;
		if(bStart<0 || bStop>e->eInfo.nBands)
			die("Index range [%d,%d) of participating bands incompatible with available bands [0,%d).\n", bStart, bStop, e->eInfo.nBands);
	}
	//Initialize spin selection:
	iSpinArr.clear();
	if(e->eInfo.spinType == SpinZ)
	{	if(spinMode==SpinAll || spinMode==SpinUp) iSpinArr.push_back(0);
		if(spinMode==SpinAll || spinMode==SpinDn) iSpinArr.push_back(1);
	}
	else
	{	assert(spinMode==SpinAll);
		iSpinArr.push_back(0);
	}
	//Initialize trial orbital centers:
	for(TrialOrbital& t: trialOrbitals)
	{	vector3<> xSum; double wSum = 0.;
		for(const Wannier::AtomicOrbital& ao: t)
		{	const vector3<>& x = (ao.sp<0 ? ao.r : e->iInfo.species[ao.sp]->atpos[ao.atom]);
			double weight = std::pow(ao.coeff, 2);
			xSum += weight * x;
			wSum += weight;
		}
		t.xCenter = (1./wSum) * xSum;
	}
	//Initialize minimizer:
	switch(localizationMeasure)
	{	case LM_FiniteDifference: wmin = std::make_shared<WannierMinimizerFD>(*e, *this); break;
		case LM_RealSpace: wmin = std::make_shared<WannierMinimizerRS>(*e, *this); break;
	}
	wmin->initTransformDependent();
	Citations::add("Maximally-localized Wannier functions",
		"N. Marzari and D. Vanderbilt, Phys. Rev. B 56, 12847 (1997)");
	
	//Initialize defects:
	for(DefectSupercell& ds: defects)
		ds.initialize(this);
}

void Wannier::saveMLWF()
{	wmin->saveMLWF();
}

string Wannier::getFilename(FilenameType fnType, string varName, int* spin) const
{	string fname;
	switch(fnType)
	{	case FilenameInit: fname = initFilename; break;
		case FilenameDump: fname = dumpFilename; break;
	}
	string spinSuffix;
	if(spin && e->eInfo.spinType==SpinZ)
		spinSuffix = (*spin)==0 ? "Up" : "Dn";
	fname.replace(fname.find("$VAR"), 4, varName+spinSuffix);
	return fname;
}

void WannierEverything::setup()
{	Everything::setup(); //base class setup
	if(!coulombParams.supercell) updateSupercell(true); //force supercell generation
	
	//Extra initialization needed for momentum, R*P and defect matrix element output:
	if(wannier.saveMomenta or wannier.saveRP or wannier.defects.size())
	{	if(eInfo.hasU)
		{	//Calculate U_rho needed for the DFT+U correction to the [r,H] momentum matrix elements:
			iInfo.rhoAtom_initZero(eVars.rhoAtom);
			iInfo.rhoAtom_initZero(eVars.U_rhoAtom);
			iInfo.rhoAtom_calc(eVars.F, eVars.C, eVars.rhoAtom);
			iInfo.rhoAtom_computeU(eVars.rhoAtom, eVars.U_rhoAtom);
		}
		bool hasUltrasoft = false;
		for(const auto& sp: iInfo.species) if(sp->isUltrasoft()) hasUltrasoft = true;
		if(hasUltrasoft or wannier.saveRP or wannier.defects.size())
		{	//Compute Vscloc needed for defects and R*P, and then call augmentDensity*Grad
			//needed for augmentation contribution to  the [r,H] momentum matrix elements:
			eVars.n = eVars.calcDensity();
			if(exCorr.needsKEdensity()) eVars.tau = eVars.KEdensity();
			eVars.EdensityAndVscloc(ener);
			iInfo.augmentDensityGridGrad(eVars.Vscloc);
		}
	}

	wannier.setup(*this);
}

