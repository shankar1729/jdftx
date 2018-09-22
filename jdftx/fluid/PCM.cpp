/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Deniz Gunceler

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

#include <fluid/PCM.h>
#include <fluid/PCM_internal.h>
#include <electronic/Everything.h>
#include <electronic/VanDerWaals.h>
#include <electronic/SpeciesInfo_internal.h>
#include <core/SphericalHarmonics.h>
#include <core/VectorField.h>
#include <core/ScalarFieldIO.h>
#include <core/Units.h>

inline double wExpand_calc(double G, double R)
{	return (2./3)*(bessel_jl(0, G*R) + bessel_jl(2, G*R)); //corresponds to theta(R-r)/(2*pi*R^3)
}

inline double wCavity_calc(double G, double d)
{	return bessel_jl(0, G*d); //corresponds to delta(d-r)
}

inline double wCavity_d_calc(double G, double d)
{	return -G*bessel_jl(1,G*d); //derivative w.r.t. d of wCavity_calc
}


//Spherically-averaged structure factor
double Sf_calc(double G, const std::vector<double>* rArr)
{	double Sf = 0.;
	for(double r: *rArr) Sf += bessel_jl(0, G*r);
	return Sf;
}


PCM::PCM(const Everything& e, const FluidSolverParams& fsp): FluidSolver(e,fsp)
{
	if(fsp.solvents.size() < 1) die("PCMs require exactly one solvent component - none specified.\n");
	if(fsp.solvents.size() > 1) die("PCMs require exactly one solvent component - more than one specified.\n");
	const auto& solvent = fsp.solvents[0];
	const double dG = 0.02;
	shape.resize(1); //single component shape by default (overridden for soft-sphere model below)
	
	//Shared information for most models:
	if(!isPCM_SCCS(fsp.pcmVariant) && (fsp.pcmVariant!=PCM_SoftSphere) && (fsp.pcmVariant!=PCM_FixedCavity))
		logPrintf("   Cavity determined by nc: %lg and sigma: %lg\n", fsp.nc, fsp.sigma);

	//Add citations, initialization and print information unique to the model:
	bool atomSpheres = false; //cavity based on atomic spheres
	bool useEta = false; //uses eta to expand electrostatic cavity
	bool cavitationNL = false; //nonlocal cavitation model
	bool dispNLatomic = false; //nonlocal dispersion model using explicit atomic sites
	bool dispNLunified = false; //nonlocal dispersion model with unified atom (single sqrtC6eff)
	double sigmaVdw = 0.; //gaussian width for unified dispersion model
	switch(fsp.pcmVariant)
	{	case PCM_SaLSA:
		{	Citations::add("Spherically-averaged liquid susceptibility ansatz (SaLSA) nonlocal fluid model",
					"R. Sundararaman, K.A. Schwarz, K. Letchworth-Weaver, and T.A. Arias, J. Chem. Phys. 142, 054102 (2015)");
			cavitationNL = true;
			dispNLatomic = true;
			break;
		}
		case PCM_CANDLE:
		{	Citations::add("Charge-asymmetric nonlocally-determined local-electric (CANDLE) solvation model",
					"R. Sundararaman and W.A. Goddard III, J. Chem. Phys. 142, 064107 (2015)");
			useEta = true;
			cavitationNL = true;
			dispNLunified = true;
			//Compute the gaussian width parameter from Rvdw:
			sigmaVdw = 1.;
			for(int iter=0; iter<50; iter++) //-- solve (Ztot wCavity * Ztot wCavity)(2 Rvdw) = nc by fixed-point (Picard) iteration
			{	double sigmaVdwNew = (2*solvent->Rvdw) / sqrt(-4. * log(fsp.nc * pow(2*sqrt(M_PI)*sigmaVdw, 3) / pow(fsp.Ztot,2)));
				if(fabs(sigmaVdwNew/sigmaVdw - 1.) < 1e-12) break;
				sigmaVdw = sigmaVdwNew;
			}
			logPrintf("   Nonlocal vdW cavity from gaussian model electron density with norm = %lg and sigma = %lg bohr\n", fsp.Ztot, sigmaVdw);
			logPrintf("   Charge asymmetry in cavity with sensitivity pCavity = %lg e-bohr/Eh\n", fsp.pCavity);
			break;
		}
		case PCM_SGA13:
		{	Citations::add("Linear/nonlinear dielectric/ionic fluid model with weighted-density cavitation and dispersion",
					"R. Sundararaman, D. Gunceler, and T.A. Arias, J. Chem. Phys. 141, 134105 (2014)");
			cavitationNL = true;
			dispNLatomic = true;
			Rex[0] = solvent->Rvdw - solvent->Res;
			Rex[1] = solvent->Rvdw;
			logPrintf("   Electrostatic cavity expanded by Rvdw-Res: %lg bohr, and cavitation/dispersion cavity by Rvdw: %lg bohr.\n", Rex[0], Rex[1]);
			//Initialize cavity expansion weight functions:
			for(int i=0; i<2; i++)
				wExpand[i].init(0, dG, e.gInfo.GmaxGrid, wExpand_calc, Rex[i]);
			break;
		}
		case PCM_GLSSA13:
		{	Citations::add("Linear/nonlinear dielectric/ionic fluid model with effective cavity tension",
				"D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, Modelling Simul. Mater. Sci. Eng. 21, 074005 (2013)");
			logPrintf("   Effective cavity tension: %lg Eh/bohr^2 to account for cavitation and dispersion.\n", fsp.cavityTension);
			break;
		}
		case PCM_LA12:
		{	if(k2factor)
				Citations::add("Linear dielectric fluid model with ionic screening",
					"K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)");
			else
				Citations::add("Linear dielectric fluid model",
					"S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)");
			logPrintf("   No cavitation or dispersion model.\n");
			break;
		}
		case PCM_SoftSphere:
		{	Citations::add("Soft-sphere continuum solvation model", "G. Fisicaro et al., J. Chem. Theory Comput. 13, 3829 (2017)");
			atomSpheres = true;
			logPrintf("   Effective cavity tension: %lg Eh/bohr^2 and pressure: %lg Eh/bohr^3 to account for cavitation and dispersion.\n", fsp.cavityTension, fsp.cavityPressure);
			logPrintf("   Cavity determined with scale factor: %lg  width: %lg  atomic radii:\n", fsp.cavityScale, fsp.sigma);
			break;
		}
		case PCM_FixedCavity:
		{	if(!fsp.cavityFile.size())
				die("Fixed-cavity solvation model requires cavityFile to be specified in pcm-params.\n\n");
			logPrintf("   Reading cavity from file '%s' ...", fsp.cavityFile.c_str()); logFlush();
			nullToZero(shape[0], gInfo);
			loadRawBinary(shape[0], fsp.cavityFile.c_str());
			fixedCavityMasked = false; //if any masks, they will be applied on first updateCavity()
			logPrintf("done.\n");
			logPrintf("   No cavitation or dispersion model.\n");
			break;
		}
		
		case_PCM_SCCS_any:
		{	Citations::add("SCCS linear dielectric fluid model", "O. Andreussi et al., J. Chem. Phys. 136, 064102 (2012)");
			if(fsp.pcmVariant==PCM_SCCS_cation || fsp.pcmVariant==PCM_SCCS_anion)
				Citations::add("SCCS parametrization for charged systems", "C. Dupont et al., J. Chem. Phys. 139, 214110 (2013)");
			logPrintf("   Cavity determined by rhoMin: %lg  rhoMax: %lg\n", fsp.rhoMin, fsp.rhoMax);
			logPrintf("   Effective cavity tension: %lg Eh/bohr^2 and pressure: %lg Eh/bohr^3 to account for cavitation and dispersion.\n",
				fsp.cavityTension, fsp.cavityPressure);
			break;
		}
	}
	
	//Initialization for atomic spheres:
	if(atomSpheres)
	{	Rall.clear();
		for(const auto& sp: e.iInfo.species)
		{	assert(sp->atomicNumber);
			double R = fsp.getAtomicRadius(*sp);
			logPrintf("\t%2s:  %.3f bohr%s\n", sp->name.c_str(), R, sp->atomicRadiusOverride ? " (Manually overriden)" : "");
			R *= fsp.cavityScale; //apply cavity scale factor
			Rall.insert(Rall.end(), sp->atpos.size(), R); //store net radius used for cavity determination
		}
		double Rmax = *std::max_element(Rall.begin(), Rall.end());
		if(fsp.ionSpacing)
		{	shape.resize(2); //need separate ionic cavity
			logPrintf("   Ionic cavity large than dielectric cavity by %lg bohr\n", fsp.ionSpacing);
			RallIonic.clear();
			for(const double R: Rall)
				RallIonic.push_back(R + fsp.ionSpacing);
			Rmax += fsp.ionSpacing;
		}
		//Initialize lattice repetitions:
		vector3<bool> isTruncated = e.coulombParams.isTruncated();
		latticeReps.clear();
		//--- circum radius of [-0.5,0.5) unit cell
		double circumR = 0.;
		vector3<int> iMax;
		for(int iDir=0; iDir<3; iDir++)
			iMax[iDir] = isTruncated[iDir] ? 0 : 1;
		vector3<int> rep;
		for(rep[0]=-iMax[0]; rep[0]<=+iMax[0]; rep[0]++)
		for(rep[1]=-iMax[1]; rep[1]<=+iMax[1]; rep[1]++)
		for(rep[2]=-iMax[2]; rep[2]<=+iMax[2]; rep[2]++)
			circumR = std::max(circumR, sqrt(gInfo.RTR.metric_length_squared(0.5*rep)));
		//--- lattice vectors in reach of this unit cell
		double inRsup = circumR + Rmax + 5.*fsp.sigma; //in-radius of supercell that contains all such points
		for(int iDir=0; iDir<3; iDir++)
			iMax[iDir] = isTruncated[iDir] ? 0 : int(ceil(inRsup * gInfo.invR.row(iDir).length()));
		for(rep[0]=-iMax[0]; rep[0]<=+iMax[0]; rep[0]++)
		for(rep[1]=-iMax[1]; rep[1]<=+iMax[1]; rep[1]++)
		for(rep[2]=-iMax[2]; rep[2]<=+iMax[2]; rep[2]++)
			if(sqrt(gInfo.RTR.metric_length_squared(rep)) <= inRsup)
				latticeReps.push_back(rep);
	}
	
	//Electrostatic cavity expansion:
	if(useEta)
	{	logPrintf("   Electrostatic cavity expanded by eta = %lg bohrs\n", fsp.eta_wDiel);
		wExpand[0].init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, wCavity_calc, fsp.eta_wDiel); //dielectric cavity expansion kernel
		wExpand[1].init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, wCavity_d_calc, fsp.eta_wDiel); //derivative of above w.r.t d
	}
	
	//Nonlocal cavitation:
	if(cavitationNL)
	{	wCavity.init(0, dG, e.gInfo.GmaxGrid, wCavity_calc, 2.*solvent->Rvdw); //Initialize nonlocal cavitation weight function
		logPrintf("   Weighted density cavitation model constrained by Nbulk: %lg bohr^-3, Pvap: %lg kPa, Rvdw: %lg bohr and sigmaBulk: %lg Eh/bohr^2 at T: %lg K.\n",
			solvent->Nbulk, solvent->Pvap/KPascal, solvent->Rvdw, solvent->sigmaBulk, fsp.T/Kelvin);
	}
	
	//Nonlocal dispersion with explicit atoms:
	if(dispNLatomic)
	{	if(!solvent->molecule.sites.size())
			die("Nonlocal dispersion model requires solvent molecule geometry, which is not yet implemented for selected solvent\n");
		Sf.resize(solvent->molecule.sites.size());
		atomicNumbers.resize(solvent->molecule.sites.size());
		for(unsigned i=0; i<Sf.size(); i++)
		{	std::vector<double> r; //radial distances of solvent sites from center
			for(vector3<> pos: solvent->molecule.sites[i]->positions) r.push_back(pos.length());
			Sf[i].init(0, dG, e.gInfo.GmaxGrid, Sf_calc, &r);
			atomicNumbers[i] = solvent->molecule.sites[i]->atomicNumber;
		}
		logPrintf("   Weighted density dispersion model using vdW pair potentials with atomic C6's and scale factor s6: %lg.\n", fsp.vdwScale);
	}
	
	//Nonlocal dispersion with united atoms:
	if(dispNLunified)
	{	logPrintf("   Weighted density dispersion model using vdW pair potentials with single solvent site with sqrtC6eff: %lg SI.\n", fsp.sqrtC6eff);
		Sf.resize(1);  //simplified model: use single site rather than explicit molecule geometry
		Sf[0].init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, RadialFunctionG::gaussTilde, 1., sigmaVdw); //CANDLE also uses this for vdW cavity
		atomicNumbers.assign(1, VanDerWaals::unitParticle); //signals point-particle with unit C6 to class VanDerWaals
	}
	
	//Optional cavity masks:
	zMask.resize(2);
	for(int iMask=0; iMask<2; iMask++)
	{	double z0 = fsp.zMask0, zH = (iMask==0 ? fsp.zMaskH : fsp.zMaskIonH);
		if(!zH) continue;
		logPrintf("   %s mask centered at z = %lg with half-width %lg in lattice coordinates, with smoothness %lg bohrs.\n",
			(iMask==0 ? "Cavity" : "Ion cavity"), z0, zH, fsp.zMaskSigma);
		if(e.coulombParams.embed)
		{	z0 -= e.coulombParams.embedCenter[2]; //offset from embedding grid center
			z0 -= floor(0.5+z0); //wrap to [-0.5,0.5)
			z0 *= 0.5; //switch to embedded double-sized grid coordinates
			zH *= 0.5;
		}
		nullToZero(zMask[iMask], gInfo);
		double* zMaskData = zMask[iMask]->data();
		const vector3<int>& S = gInfo.S;
		const size_t iStart=0, iStop = gInfo.nr;
		double SzInv = 1./S[2];
		double sigmaInv = gInfo.R.column(2).length()/fsp.zMaskSigma; //1/sigma in lattice coordinates
		THREAD_rLoop(
			double zDiff = SzInv*iv[2] - z0;
			zDiff -= floor(0.5+zDiff); //wrap to [-0.5,0.5)
			zMaskData[i] = 0.5*erfc(sigmaInv*(zH - fabs(zDiff))); //smoothly goes to 0 when |zDiff| < zH
		)
	}
	nShape = shape.size(); //will be 2 when separate ionic cavity is used, 1 otherwise
}

PCM::~PCM()
{	for(int i=0; i<2; i++) wExpand[i].free();
	wCavity.free();
	for(unsigned i=0; i<Sf.size(); i++) Sf[i].free();
}

void PCM::updateCavity()
{
	bool cavityChanged = true; //keep track of whether cavity is updated in code below (usually yes, except for SS)
	
	//Cavities from expanded densities for SGA13 variant:
	if(fsp.pcmVariant == PCM_SGA13)
	{	ScalarField* shapeEx[2] = { &shape[0], &shapeVdw };
		for(int i=0; i<2; i++)
		{	ShapeFunctionSGA13::expandDensity(wExpand[i], Rex[i], nCavity, nCavityEx[i]);
			ShapeFunction::compute(nCavityEx[i], *(shapeEx[i]), fsp.nc, fsp.sigma);
		}
	}
	else if(fsp.pcmVariant == PCM_CANDLE)
	{	nCavityEx[0] = fsp.Ztot * I(Sf[0] * J(nCavity));
		ShapeFunctionCANDLE::compute(nCavityEx[0], coulomb(Sf[0]*rhoExplicitTilde), shapeVdw,
			fsp.nc, fsp.sigma, fsp.pCavity); //vdW cavity
		shape[0] = I(wExpand[0] * J(shapeVdw)); //dielectric cavity
	}
	else if(fsp.pcmVariant==PCM_SoftSphere)
	{	//Construct flat list of all atom positions:
		std::vector<vector3<>> atposAllCur;
		for(const std::vector<vector3<>>& atposSp: atpos)
			atposAllCur.insert(atposAllCur.end(), atposSp.begin(), atposSp.end());
		//Update cavity only if atomic positions have changed:
		if(atposAllCur != atposAll)
		{	atposAll = atposAllCur;
			nullToZero(shape, gInfo);
			ShapeFunctionSoftSphere::compute(atposAll, latticeReps, Rall, shape[0], fsp.sigma);
			if(fsp.ionSpacing)
				ShapeFunctionSoftSphere::compute(atposAll, latticeReps, RallIonic, shape[1], fsp.sigma);
		}
		else cavityChanged = false; //cavity not updated (so don't re-apply masks below)
	}
	else if(fsp.pcmVariant == PCM_FixedCavity)
	{	//Cavity read in at the beginning and never changed subsequently
		cavityChanged = (not fixedCavityMasked); //true on first call, so that mask applied below
		fixedCavityMasked = true; //mask will get applied this time, so need not be done again after.
	}
	else if(isPCM_SCCS(fsp.pcmVariant))
		ShapeFunctionSCCS::compute(nCavity, shape[0], fsp.rhoMin, fsp.rhoMax, epsBulk);
	else //Compute directly from nCavity (which is a density product for SaLSA):
		ShapeFunction::compute(nCavity, shape[0], fsp.nc, fsp.sigma);
	
	//Apply cavity masks (if any):
	if((zMask[0] || zMask[1]) && cavityChanged)
	{	if(nShape==1)
		{	if(zMask[1]) //ion mask specified on model with natually only one cavity
			{	shape.resize(2); //so make ionic cavity separate based on specified ion mask
				shape[1] = zMask[1] * shape[0];
			}
			if(zMask[0]) shape[0] *= zMask[0];
		}
		else //separate cavities anyway: apply masks separately
			for(int iMask=0; iMask<2; iMask++)
				if(zMask[iMask])
					shape[iMask] *= zMask[iMask];
	}
	
	//Compute and cache cavitation energy and gradients:
	const auto& solvent = fsp.solvents[0];
	switch(fsp.pcmVariant)
	{	case PCM_SaLSA:
		case PCM_CANDLE:
		case PCM_SGA13:
		{	//Select relevant shape function:
			const ScalarFieldTilde sTilde = J(fsp.pcmVariant==PCM_SaLSA ? shape[0] : shapeVdw);
			ScalarFieldTilde A_sTilde;
			//Cavitation:
			const double nlT = solvent->Nbulk * fsp.T;
			const double Gamma = log(nlT/solvent->Pvap) - 1.;
			const double Cp = 15. * (solvent->sigmaBulk/(2*solvent->Rvdw * nlT) - (1+Gamma)/6);
			const double coeff2 = 1. + Cp - 2.*Gamma;
			const double coeff3 = Gamma - 1. -2.*Cp;
			ScalarField sbar = I(wCavity*sTilde);
			Adiel["Cavitation"] = nlT * integral(sbar*(Gamma + sbar*(coeff2 + sbar*(coeff3 + sbar*Cp))));
			A_sTilde += wCavity*Idag(nlT * (Gamma + sbar*(2.*coeff2 + sbar*(3.*coeff3 + sbar*(4.*Cp)))));
			//Dispersion:
			ScalarFieldTildeArray Ntilde(Sf.size()), A_Ntilde(Sf.size()); //effective nuclear densities in spherical-averaged ansatz
			for(unsigned i=0; i<Sf.size(); i++)
				Ntilde[i] = solvent->Nbulk * (Sf[i] * sTilde);
			const double vdwScaleEff = (fsp.pcmVariant==PCM_CANDLE) ? fsp.sqrtC6eff : fsp.vdwScale;
			Adiel["Dispersion"] = e.vanDerWaals->energyAndGrad(atpos, Ntilde, atomicNumbers, vdwScaleEff, &A_Ntilde);
			A_vdwScale = Adiel["Dispersion"]/vdwScaleEff;
			for(unsigned i=0; i<Sf.size(); i++)
				if(A_Ntilde[i])
					A_sTilde += solvent->Nbulk * (Sf[i] * A_Ntilde[i]);
			//Propagate gradients to appropriate shape function:
			(fsp.pcmVariant==PCM_SaLSA ? Acavity_shape : Acavity_shapeVdw) = Jdag(A_sTilde);
			break;
		}
		case PCM_GLSSA13:
		case PCM_SoftSphere:
		{	VectorField Dshape = gradient(shape[0]);
			ScalarField surfaceDensity = sqrt(lengthSquared(Dshape));
			ScalarField invSurfaceDensity = inv(surfaceDensity);
			A_tension = integral(surfaceDensity);
			Adiel["CavityTension"] = A_tension * fsp.cavityTension;
			Acavity_shape = (-fsp.cavityTension)*divergence(Dshape*invSurfaceDensity);
			if(fsp.cavityPressure)
			{	Adiel["CavityPressure"] = fsp.cavityPressure * (gInfo.detR - integral(shape[0]));
				Acavity_shape = Acavity_shape - fsp.cavityPressure;
			}
			break;
		}
		case PCM_FixedCavity:
		case PCM_LA12:
			break; //no contribution
		case_PCM_SCCS_any:
		{	//Volume contribution:
			Adiel["CavityPressure"] = fsp.cavityPressure * (gInfo.detR - integral(shape[0]));
			//Surface contribution:
			ScalarField shapePlus, shapeMinus;
			ShapeFunctionSCCS::compute(nCavity+(0.5*fsp.rhoDelta), shapePlus, fsp.rhoMin, fsp.rhoMax, epsBulk);
			ShapeFunctionSCCS::compute(nCavity-(0.5*fsp.rhoDelta), shapeMinus, fsp.rhoMin, fsp.rhoMax, epsBulk);
			ScalarField DnLength = sqrt(lengthSquared(gradient(nCavity)));
			Adiel["CavityTension"] = (fsp.cavityTension/fsp.rhoDelta) * integral(DnLength * (shapeMinus - shapePlus));
			break;
		}
	}
}

void PCM::propagateCavityGradients(const ScalarFieldArray& A_shape, ScalarField& A_nCavity, ScalarFieldTilde& A_rhoExplicitTilde, IonicGradient* forces) const
{
	if(forces) forces->init(e.iInfo); //zero and initialize forces if needed
	
	//Account for cavity masks, if any
	if(zMask[0] || zMask[1])
	{	if(nShape==1)
		{	if(zMask[0]) (ScalarField&)A_shape[0] *= zMask[0];
			if(zMask[1]) (ScalarField&)A_shape[0] += zMask[1] * A_shape[1]; //shape[1] was created from shape[0] above
		}
		else //separate cavities anyway: apply masks separately
			for(int iMask=0; iMask<2; iMask++)
				if(zMask[iMask])
					(ScalarField&)A_shape[iMask] *= zMask[iMask];
	}
	
	if(fsp.pcmVariant == PCM_SGA13)
	{	//Propagate gradient w.r.t expanded cavities to nCavity:
		((PCM*)this)->A_nc = 0;
		const ScalarField* A_shapeEx[2] = { &A_shape[0], &Acavity_shapeVdw };
		for(int i=0; i<2; i++)
		{	//First compute derivative w.r.t expanded electron density:
			ScalarField A_nCavityEx;
			ShapeFunction::propagateGradient(nCavityEx[i], *(A_shapeEx[i]), A_nCavityEx, fsp.nc, fsp.sigma);
			((PCM*)this)->A_nc += (-1./fsp.nc) * integral(A_nCavityEx*nCavityEx[i]);
			//then propagate to original electron density:
			ScalarField nCavityExUnused; //unused return value below
			ShapeFunctionSGA13::expandDensity(wExpand[i], Rex[i], nCavity, nCavityExUnused, &A_nCavityEx, &A_nCavity);
		}
	}
	else if(fsp.pcmVariant == PCM_CANDLE)
	{	ScalarField A_nCavityEx; ScalarFieldTilde A_phiExt; double A_pCavity=0.;
		ShapeFunctionCANDLE::propagateGradient(nCavityEx[0], coulomb(Sf[0]*rhoExplicitTilde), I(wExpand[0]*J(A_shape[0])) + Acavity_shapeVdw,
			A_nCavityEx, A_phiExt, A_pCavity, fsp.nc, fsp.sigma, fsp.pCavity);
		A_nCavity += fsp.Ztot * I(Sf[0] * J(A_nCavityEx));
		((PCM*)this)->A_rhoNonES = coulomb(Sf[0]*A_phiExt);
		A_rhoExplicitTilde += A_rhoNonES;
		((PCM*)this)->A_nc = (-1./fsp.nc) * integral(A_nCavityEx*nCavityEx[0]);
		((PCM*)this)->A_eta_wDiel = integral(A_shape[0] * I(wExpand[1]*J(shapeVdw)));
		((PCM*)this)->A_pCavity = A_pCavity;
	}
	else if(fsp.pcmVariant == PCM_SoftSphere)
	{	nullToZero(A_nCavity, gInfo); //no electronic contributions
		if(forces)
		{	std::vector<vector3<>> A_atposAll; //gradient w.r.t atomic positions
			std::vector<double> A_Rall; //gradient w.r.t atomic radii
			ShapeFunctionSoftSphere::propagateGradient(atposAll, latticeReps, Rall, shape[0], A_shape[0]+Acavity_shape, A_atposAll, A_Rall, fsp.sigma);
			if(fsp.ionSpacing)
				ShapeFunctionSoftSphere::propagateGradient(atposAll, latticeReps, RallIonic, shape[1], A_shape[1], A_atposAll, A_Rall, fsp.sigma);
			auto A_atposAllPtr = A_atposAll.begin();
			for(unsigned iSp=0; iSp<atpos.size(); iSp++)
				for(unsigned iAtom=0; iAtom<atpos[iSp].size(); iAtom++)
					(*forces)[iSp][iAtom] -= *(A_atposAllPtr++); //negative gradient
			((PCM*)this)->A_cavityScale = cblas_ddot(Rall.size(), Rall.data(),1, A_Rall.data(),1) / fsp.cavityScale; //gradient w.r.t scale factor used for fits
		}
	}
	else if(fsp.pcmVariant == PCM_FixedCavity)
	{	nullToZero(A_nCavity, gInfo); //No electronic or force contributions
	}
	else if(isPCM_SCCS(fsp.pcmVariant))
	{	//Electrostatic and volumetric combinations via shape:
		ShapeFunctionSCCS::propagateGradient(nCavity, A_shape[0] - fsp.cavityPressure, A_nCavity, fsp.rhoMin, fsp.rhoMax, epsBulk);
		//Add surface contributions:
		ScalarField shapePlus, shapeMinus;
		ShapeFunctionSCCS::compute(nCavity+(0.5*fsp.rhoDelta), shapePlus, fsp.rhoMin, fsp.rhoMax, epsBulk);
		ShapeFunctionSCCS::compute(nCavity-(0.5*fsp.rhoDelta), shapeMinus, fsp.rhoMin, fsp.rhoMax, epsBulk);
		VectorField Dn = gradient(nCavity);
		ScalarField DnLength = sqrt(lengthSquared(Dn));
		ScalarField A_shapeMinus = (fsp.cavityTension/fsp.rhoDelta) * DnLength;
		ScalarField A_DnLength = (fsp.cavityTension/fsp.rhoDelta) * (shapeMinus - shapePlus);
		A_nCavity -= divergence(Dn * (inv(DnLength) * A_DnLength));
		ShapeFunctionSCCS::propagateGradient(nCavity+(0.5*fsp.rhoDelta), -A_shapeMinus, A_nCavity, fsp.rhoMin, fsp.rhoMax, epsBulk);
		ShapeFunctionSCCS::propagateGradient(nCavity-(0.5*fsp.rhoDelta),  A_shapeMinus, A_nCavity, fsp.rhoMin, fsp.rhoMax, epsBulk);
	}
	else //All gradients are w.r.t the same shape function - propagate them to nCavity (which is defined as a density product for SaLSA)
	{	ShapeFunction::propagateGradient(nCavity, A_shape[0] + Acavity_shape, A_nCavity, fsp.nc, fsp.sigma);
		((PCM*)this)->A_nc = (-1./fsp.nc) * integral(A_nCavity*nCavity);
	}
}

void PCM::accumExtraForces(IonicGradient* forces, const ScalarFieldTilde& A_nCavityTilde) const
{	if(forces)
	{
		//VDW contribution:
		switch(fsp.pcmVariant)
		{	case PCM_SaLSA:
			case PCM_CANDLE:
			case PCM_SGA13:
			{	const auto& solvent = fsp.solvents[0];
				const ScalarFieldTilde sTilde = J(fsp.pcmVariant==PCM_SaLSA ? shape[0] : shapeVdw);
				ScalarFieldTildeArray Ntilde(Sf.size());
				for(unsigned i=0; i<Sf.size(); i++)
					Ntilde[i] = solvent->Nbulk * (Sf[i] * sTilde);
				const double vdwScaleEff = (fsp.pcmVariant==PCM_CANDLE) ? fsp.sqrtC6eff : fsp.vdwScale;
				e.vanDerWaals->energyAndGrad(atpos, Ntilde, atomicNumbers, vdwScaleEff, 0, forces);
				break;
			}
			default: break; //no VDW forces
		}
		//Full core contribution:
		switch(fsp.pcmVariant)
		{	case PCM_SaLSA:
			case PCM_CANDLE:
			{	VectorFieldTilde gradAtpos; nullToZero(gradAtpos, gInfo);
				for(unsigned iSp=0; iSp<atpos.size(); iSp++)
					for(unsigned iAtom=0; iAtom<atpos[iSp].size(); iAtom++)
					{	callPref(gradSGtoAtpos)(gInfo.S, atpos[iSp][iAtom], A_nCavityTilde->dataPref(), gradAtpos.dataPref());
						for(int k=0; k<3; k++)
							(*forces)[iSp][iAtom][k] -= e.iInfo.species[iSp]->ZfullCore * sum(gradAtpos[k]); //negative gradient
					}
				break;
			}
			default: break; //no full-core forces
		}
	}
}

ScalarFieldTilde PCM::getFullCore() const
{	switch(fsp.pcmVariant)
	{	case PCM_SaLSA:
		case PCM_CANDLE:
		{	ScalarFieldTilde nFullCore, SG(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
			for(unsigned iSp=0; iSp<atpos.size(); iSp++)
			{	ManagedArray<vector3<>> atposManaged(atpos[iSp]);
				//Compute structure factor and accumulate contribution from this species:
				callPref(getSG)(gInfo.S, atposManaged.nData(), atposManaged.dataPref(), 1./gInfo.detR, SG->dataPref());
				nFullCore += e.iInfo.species[iSp]->ZfullCore * SG;
			}
			return nFullCore;
		}
		default: return 0; //no full-core
	}
}


void PCM::dumpDensities(const char* filenamePattern) const
{	string filename;
	FLUID_DUMP(shape[0], "Shape");
    if(shape.size() > 1)
	{	FLUID_DUMP(shape[1], "ShapeIonic");
	}
	if(fsp.pcmVariant==PCM_SGA13 || fsp.pcmVariant==PCM_CANDLE)
	{	FLUID_DUMP(shapeVdw, "ShapeVdw");
	}
}

void PCM::dumpDebug(const char* filenamePattern) const
{	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Debug");
	logPrintf("Dumping '%s' ... ", filename.c_str());  logFlush();
	FILE* fp = mpiWorld->isHead() ? fopen(filename.c_str(), "w") : nullLog;
	if(!fp) die("Error opening %s for writing.\n", filename.c_str());

	fprintf(fp, "Dielectric cavity volume = %f\n", integral(1.-shape[0]));
	fprintf(fp, "Dielectric cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shape[0])))));
	if(shape.size() > 1)
	{	fprintf(fp, "Ionic cavity volume = %f\n", integral(1.-shape[1]));
		fprintf(fp, "Ionic cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shape[1])))));
	}
	if(fsp.pcmVariant==PCM_SGA13 || fsp.pcmVariant==PCM_CANDLE)
	{	fprintf(fp, "VDW cavity volume = %f\n", integral(1.-shapeVdw));
		fprintf(fp, "VDW cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shapeVdw)))));
	}
	
	fprintf(fp, "\nComponents of Adiel:\n");
	Adiel.print(fp, true, "   %13s = %25.16lf\n");	
	
	fprintf(fp, "\n\nGradients wrt fit parameters:\n");
	if((!isPCM_SCCS(fsp.pcmVariant)) && (fsp.pcmVariant!=PCM_SoftSphere))
		fprintf(fp, "   E_nc = %.15lg\n", A_nc);
	switch(fsp.pcmVariant)
	{	case PCM_SaLSA:
		case PCM_SGA13:
			fprintf(fp, "   E_vdwScale = %.15lg\n", A_vdwScale);
			break;
		case PCM_CANDLE:
			fprintf(fp, "   E_sqrtC6eff = %.15lg\n", A_vdwScale);
			fprintf(fp, "   E_eta_wDiel = %.15lg\n", A_eta_wDiel);
			fprintf(fp, "   E_pCavity = %.15lg\n", A_pCavity);
			break;
		case PCM_GLSSA13:
			fprintf(fp, "   E_tension = %.15lg\n", A_tension);
			break;
		case PCM_SoftSphere:
			fprintf(fp, "   E_cavityScale = %.15lg\n", A_cavityScale);
			fprintf(fp, "   E_tension = %.15lg\n", A_tension);
		case PCM_FixedCavity:
		case PCM_LA12:
		case_PCM_SCCS_any:
			break;
	}
	printDebug(fp);

	if(mpiWorld->isHead()) fclose(fp);
	logPrintf("done\n"); logFlush();
	
	{ //scope for overriding filename
		char filename[256];	ostringstream oss;
		oss << "Nspherical";
		sprintf(filename, filenamePattern, oss.str().c_str());
		logPrintf("Dumping '%s' ... ", filename); logFlush();
		saveSphericalized(&shape[0], 1, filename);
		logPrintf("done\n"); logFlush();
	}
	
	if(fsp.pcmVariant==PCM_SGA13 || fsp.pcmVariant==PCM_CANDLE )
	{
		char filename[256];	ostringstream oss;
		oss << "NvdWspherical";
		sprintf(filename, filenamePattern, oss.str().c_str());
		logPrintf("Dumping '%s' ... ", filename); logFlush();
		saveSphericalized(&shapeVdw,1, filename);
		logPrintf("done\n"); logFlush();
	}
}
