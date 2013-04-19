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

#include <electronic/PCM.h>
#include <electronic/PCM_internal.h>
#include <electronic/Everything.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/VanDerWaals.h>
#include <core/DataMultiplet.h>
#include <core/DataIO.h>

void wExpand_calc(int i, double G2, double R, double* w)
{	double GR = sqrt(G2) * R;
	w[i] = (2./3)*(bessel_jl(0, GR) + bessel_jl(2, GR)); //corresponds to theta(R-r)/(2*pi*R^3)
}

void wCavity_calc(int i, double G2, double d, double* w)
{	w[i] = bessel_jl(0, sqrt(G2) * d); //corresponds to delta(d-r)
}

//Spherically-averaged structure factor
void Sf_calc(int i, double G2, const std::vector<double>* rArr, double* Sf)
{	double G = sqrt(G2);
	Sf[i] = 0.;
	for(double r: *rArr)
		Sf[i] += bessel_jl(0, G*r);
}


PCM::PCM(const Everything& e, const FluidSolverParams& fsp): FluidSolver(e), fsp(fsp)
{
	k2factor = (8*M_PI/fsp.T) * fsp.ionicConcentration * pow(fsp.ionicZelectrolyte,2);

	//Print common info and add relevant citations:
	logPrintf("   Cavity determined by nc: %lg and sigma: %lg\n", fsp.nc, fsp.sigma);
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
			Citations::add("Linear/nonlinear dielectric/ionic fluid model with weighted-density cavitation and dispersion",
				"R. Sundararaman, D. Gunceler, and T.A. Arias, (under preparation)");
			Rex[0] = fsp.Rvdw -fsp.Res;
			Rex[1] = fsp.Rvdw;
			logPrintf("   Electrostatic cavity expanded by Rvdw-Res: %lg bohr, and cavitation/dispersion cavity by Rvdw: %lg bohr.\n", Rex[0], Rex[1]);
			logPrintf("   Weighted density cavitation model constrained by Nbulk: %lg bohr^-3, Pvap: %lg kPa and sigmaBulk: %lg Eh/bohr^2 at T: %lg K.\n", fsp.Nbulk, fsp.Pvap/KPascal, fsp.sigmaBulk, fsp.T/Kelvin);
			logPrintf("   Weighted density dispersion model using vdW pair potentials.\n");
			//Initialize cavity expansion weight functions:
			for(int i=0; i<2; i++)
			{	wExpand[i] = std::make_shared<RealKernel>(e.gInfo);
				applyFuncGsq(e.gInfo, wExpand_calc, Rex[i], wExpand[i]->data);
				wExpand[i]->set();
			}
			//Initialize nonlocal cavitation weight function:
			wCavity = std::make_shared<RealKernel>(e.gInfo);
			applyFuncGsq(e.gInfo, wCavity_calc, 2.*fsp.Rvdw, wCavity->data);
			wCavity->set();
			//Initialize structure factors for dispersion:
			if(!fsp.pcmSite.size()) die("Nonlocal dispersion model requires solvent molecule geometry, which is not yet implemented for selected solvent\n");
			Sf.resize(fsp.pcmSite.size());
			for(unsigned i=0; i<Sf.size(); i++)
			{	std::vector<double> r; //radial distances of solvent sites from center
				for(vector3<> pos: fsp.pcmSite[i].pos) r.push_back(pos.length());
				Sf[i] = std::make_shared<RealKernel>(e.gInfo);
				applyFuncGsq(e.gInfo, Sf_calc, &r, Sf[i]->data);
				Sf[i]->set();
			}
			vdwForces = std::make_shared<IonicGradient>();
			break;
		case PCM_GLSSA13:
			Citations::add("Linear/nonlinear dielectric/ionic fluid model with effective cavity tension",
				"D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189");
			logPrintf("   Effective cavity tension: %lg Eh/bohr^2 to account for cavitation and dispersion.\n", fsp.cavityTension);
			break;
		case PCM_LA12:
		case PCM_PRA05:
			if(fsp.ionicConcentration)
				Citations::add("Linear dielectric fluid model with ionic screening",
					"K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)");
			else
				Citations::add("Linear dielectric fluid model",
					"S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)");
			logPrintf("   No cavitation model.\n");
			break;
	}
}

void PCM::updateCavity()
{
	//Cavities from expanded densities for SGA13 variant:
	if(fsp.pcmVariant == PCM_SGA13)
	{	DataRptr* shapeEx[2] = { &shape, &shapeVdw };
		for(int i=0; i<2; i++)
		{	ShapeFunction::expandDensity(*(wExpand[i]), Rex[i], nCavity, nCavityEx[i]);
			ShapeFunction::compute(nCavityEx[i], *(shapeEx[i]), fsp.nc, fsp.sigma);
		}
	}
	else //Cavity from actual electron density:
		ShapeFunction::compute(nCavity, shape, fsp.nc, fsp.sigma);
	
	//Compute and cache cavitation energy and gradients:
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
		{	//Select relevant shape function:
			const DataGptr sTilde = J(fsp.pcmVariant==PCM_SGA13 ? shapeVdw : shape);
			DataGptr A_sTilde;
			//Cavitation:
			const double nlT = fsp.Nbulk * fsp.T;
			const double Gamma = log(nlT/fsp.Pvap) - 1.;
			const double Cp = 15. * (fsp.sigmaBulk/(2*fsp.Rvdw * nlT) - (1+Gamma)/6);
			const double coeff2 = 1. + Cp - 2.*Gamma;
			const double coeff3 = Gamma - 1. -2.*Cp;
			DataRptr sbar = I((*wCavity)*sTilde);
			Adiel["Cavitation"] = nlT * integral(sbar*(Gamma + sbar*(coeff2 + sbar*(coeff3 + sbar*Cp))));
			A_sTilde += (*wCavity)*Idag(nlT * (Gamma + sbar*(2.*coeff2 + sbar*(3.*coeff3 + sbar*(4.*Cp)))));
			//Dispersion:
			DataGptrCollection Ntilde(Sf.size()), A_Ntilde(Sf.size()); //effective nuclear densities in spherical-averaged ansatz
			std::vector<int> atomicNumbers(Sf.size());
			for(unsigned i=0; i<Sf.size(); i++)
			{	Ntilde[i] = fsp.Nbulk * ((*Sf[i]) * sTilde);
				atomicNumbers[i] = fsp.pcmSite[i].atomicNumber;
			}
			vdwForces->init(e.iInfo);
			Adiel["Dispersion"] = e.vanDerWaals->energyAndGrad(Ntilde, atomicNumbers, fsp.VDWCouplingScale, &A_Ntilde, &(*vdwForces));
			A_vdwScale = Adiel["Dispersion"]/fsp.VDWCouplingScale;
			for(unsigned i=0; i<Sf.size(); i++)
				A_sTilde += fsp.Nbulk * ((*Sf[i]) * A_Ntilde[i]);
			//Propagate gradients to appropriate shape function:
			(fsp.pcmVariant==PCM_SGA13 ? Acavity_shapeVdw : Acavity_shape) = Jdag(A_sTilde);
			break;
		}
		case PCM_GLSSA13:
		{	DataRptrVec Dshape = gradient(shape);
			DataRptr surfaceDensity = sqrt(lengthSquared(Dshape));
			DataRptr invSurfaceDensity = inv(surfaceDensity);
			A_tension = integral(surfaceDensity);
			Adiel["CavityTension"] = A_tension * fsp.cavityTension;
			Acavity_shape = (-fsp.cavityTension)*divergence(Dshape*invSurfaceDensity);
			break;
		}
		case PCM_LA12:
		case PCM_PRA05:
		default:
			break; //no contribution
	}
}

void PCM::propagateCavityGradients(const DataRptr& A_shape, DataRptr& A_nCavity) const
{	if(fsp.pcmVariant == PCM_SGA13)
	{	//Propagate gradient w.r.t expanded cavities to nCavity:
		A_nCavity = 0;
		((PCM*)this)->A_nc = 0;
		const DataRptr* A_shapeEx[2] = { &A_shape, &Acavity_shapeVdw };
		for(int i=0; i<2; i++)
		{	//First compute derivative w.r.t expanded electron density:
			DataRptr A_nCavityEx;
			ShapeFunction::propagateGradient(nCavityEx[i], *(A_shapeEx[i]), A_nCavityEx, fsp.nc, fsp.sigma);
			((PCM*)this)->A_nc += (-1./fsp.nc) * integral(A_nCavityEx*nCavityEx[i]);
			//then propagate to original electron density:
			DataRptr nCavityExUnused; //unused return value below
			ShapeFunction::expandDensity(*(wExpand[i]), Rex[i], nCavity, nCavityExUnused, &A_nCavityEx, &A_nCavity);
		}
	}
	else //All gradients are w.r.t the same shape function - propagate them to nCavity
	{	A_nCavity = 0;
		ShapeFunction::propagateGradient(nCavity, A_shape + Acavity_shape, A_nCavity, fsp.nc, fsp.sigma);
		((PCM*)this)->A_nc = (-1./fsp.nc) * integral(A_nCavity*nCavity);
	}
}

void PCM::dumpDensities(const char* filenamePattern) const
{	string filename;
	FLUID_DUMP(shape, "Shape");
    if(fsp.pcmVariant == PCM_SGA13)
	{	FLUID_DUMP(shapeVdw, "ShapeVdw");
	}
}

void PCM::dumpDebug(const char* filenamePattern) const
{	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Debug");
	logPrintf("Dumping '%s'... \t", filename.c_str());  logFlush();
	FILE* fp = fopen(filename.c_str(), "w");
	if(!fp) die("Error opening %s for writing.\n", filename.c_str());	

	fprintf(fp, "Cavity volume = %f\n", integral(1.-shape));
	fprintf(fp, "Cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shape)))));
	if(fsp.pcmVariant == PCM_SGA13)
	{	fprintf(fp, "Expanded cavity volume = %f\n", integral(1.-shapeVdw));
		fprintf(fp, "Expanded cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shapeVdw)))));
	}
	
	fprintf(fp, "\nComponents of Adiel:\n");
	Adiel.print(fp, true, "   %13s = %25.16lf\n");	
	
	fprintf(fp, "\n\nGradients wrt fit parameters:\n");
	fprintf(fp, "   E_nc = %f\n", A_nc);
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
			fprintf(fp, "   E_vdwScale = %f\n", A_vdwScale);
			break;
		case PCM_GLSSA13:
			fprintf(fp, "   E_t = %f\n", A_tension);
			break;
		case PCM_LA12:
		case PCM_PRA05:
			break;
	}

	printDebug(fp);
	
	fclose(fp);
	logPrintf("done\n"); logFlush();
}
