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
#include <core/DataMultiplet.h>
#include <core/DataIO.h>

void wExpand_calc(int i, double G2, double R, double* w)
{	double GR = sqrt(G2) * R;
	w[i] = (2./3)*(bessel_jl(0, GR) + bessel_jl(2, GR)); //corresponds to theta(R-r)/(2*pi*R^3)
}

void wCavity_calc(int i, double G2, double d, double* w)
{	w[i] = bessel_jl(0, sqrt(G2) * d); //corresponds to delta(d-r)
}

PCM::PCM(const Everything& e, const FluidSolverParams& fsp): FluidSolver(e), fsp(fsp), wExpand(0), wCavity(0)
{
	k2factor = (8*M_PI/fsp.T) * fsp.ionicConcentration * pow(fsp.ionicZelectrolyte,2);

	//Print common info and add relevant citations:
	logPrintf("   Cavity determined by nc: %lg and sigma: %lg\n", fsp.nc, fsp.sigma);
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
			Citations::add("Linear/nonlinear dielectric/ionic fluid model with weighted-density cavitation and dispersion",
				"R. Sundararaman, D. Gunceler, and T.A. Arias, (under preparation)");
			logPrintf("   Weighted density cavitation model constrained by Nbulk: %lg bohr^-3, Pvap: %lg kPa and sigmaBulk: %lg Eh/bohr^2 at T: %lg K.\n", fsp.Nbulk, fsp.Pvap/KPascal, fsp.sigmaBulk, fsp.T/Kelvin);
			logPrintf("   Weighted density dispersion model using vdW pair potentials.\n");
			//Initialize cavity expansion weight function:
			wExpand = new RealKernel(e.gInfo);
			applyFuncGsq(e.gInfo, wExpand_calc, fsp.Rvdw, wExpand->data);
			wExpand->set();
			//Initialize nonlocal cavitation weight function:
			wCavity = new RealKernel(e.gInfo);
			applyFuncGsq(e.gInfo, wCavity_calc, 2.*fsp.Rvdw, wCavity->data);
			wCavity->set();
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

PCM::~PCM()
{	if(wExpand) delete wExpand;
	if(wCavity) delete wCavity;
}

void PCM::updateCavity()
{
	//Compute cavity shape function
	ShapeFunction::compute(nCavity, shape, fsp.nc, fsp.sigma);
	
	//Expanded cavity for SGA13 variant:
	if(fsp.pcmVariant == PCM_SGA13)
	{	ShapeFunction::expandDensity(*wExpand, fsp.Rvdw, nCavity, nCavityEx);
		ShapeFunction::compute(nCavityEx, shapeEx, fsp.nc, fsp.sigma);
	}
	
	//Compute and cache cavitation energy and gradients:
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
		{	//Cavitation:
			const double nlT = fsp.Nbulk * fsp.T;
			const double Gamma = log(nlT/fsp.Pvap) - 1.;
			const double Cp = 15. * (fsp.sigmaBulk/(2*fsp.Rvdw * nlT) - (1+Gamma)/6);
			const double coeff2 = 1. + Cp - 2.*Gamma;
			const double coeff3 = Gamma - 1. -2.*Cp;
			DataRptr sbar = I((*wCavity)*J(shapeEx));
			Adiel["Cavitation"] = nlT * integral(sbar*(Gamma + sbar*(coeff2 + sbar*(coeff3 + sbar*Cp))));
			Acavity_shapeEx = Jdag((*wCavity)*Idag(nlT * (Gamma + sbar*(2.*coeff2 + sbar*(3.*coeff3 + sbar*(4.*Cp))))));
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
{	//Compute nCavity gradient including cached cavitation gradient:
	A_nCavity = 0;
	ShapeFunction::propagateGradient(nCavity, A_shape + Acavity_shape, A_nCavity, fsp.nc, fsp.sigma);
	((PCM*)this)->A_nc = (-1./fsp.nc) * integral(A_nCavity*nCavity);
	
	//Propagate cavitation/dispersion gradient w.r.t expanded cavity:
	if(fsp.pcmVariant == PCM_SGA13)
	{	//First compute derivative w.r.t expanded electron density:
		DataRptr A_nCavityEx;
		ShapeFunction::propagateGradient(nCavityEx, Acavity_shapeEx, A_nCavityEx, fsp.nc, fsp.sigma);
		((PCM*)this)->A_nc += (-1./fsp.nc) * integral(A_nCavityEx*nCavityEx);
		//then propagate to original electron density:
		DataRptr nCavityExUnused; //unused return value below
		ShapeFunction::expandDensity(*wExpand, fsp.Rvdw, nCavity, nCavityExUnused, &A_nCavityEx, &A_nCavity);
	}
}

void PCM::dumpDensities(const char* filenamePattern) const
{	string filename;
	FLUID_DUMP(shape, "Shape");
    if(fsp.pcmVariant == PCM_SGA13)
	{	FLUID_DUMP(shapeEx, "ShapeEx");
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
	{	fprintf(fp, "Expanded cavity volume = %f\n", integral(1.-shapeEx));
		fprintf(fp, "Expanded cavity surface area = %f\n", integral(sqrt(lengthSquared(gradient(shapeEx)))));
	}
	
	fprintf(fp, "\nComponents of Adiel:\n");
	Adiel.print(fp, true, "   %13s = %25.16lf\n");	
	
	fprintf(fp, "\n\nGradients wrt fit parameters:\n");
	fprintf(fp, "   E_nc = %f\n", A_nc);
	fprintf(fp, "   E_t = %f\n", A_tension);

	printDebug(fp);
	
	fclose(fp);
	logPrintf("done\n"); logFlush();
}
