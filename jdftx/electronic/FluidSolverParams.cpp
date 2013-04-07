/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#include <electronic/FluidSolverParams.h>
#include <electronic/symbols.h>
#include <fluid/Fex_H2O_ScalarEOS_internal.h>

FluidSolverParams::FluidSolverParams()
: verboseLog(false),
linearDielectric(false), linearScreening(false),
VDWCouplingScale(0.75)
{
}

//! Vapor pressure from the Antoine equation
//! @param T temperature (in a.u.)
//! @param A log prefactor for pressure in KPascal
//! @param B in Kelvin
//! @param C in Kelvin
inline double antoinePvap(double T, double A, double B, double C)
{	return KPascal * pow(10, A - B*Kelvin/(C*Kelvin + T));
}

//Wrapper to atomicSymbolMap
inline int getAtomicNumber(const char* symbol)
{	int atNum = 0;
	atomicSymbolMap.getEnum(symbol, atNum);
	return atNum;
}

void FluidSolverParams::setPCMparams()
{
	//Set physical parameters (in atomic units) decsribing solvent:
	switch(solventName)
	{	case H2O:
		{	epsBulk = 78.4;
			Nbulk = 4.9383e-3;
			pMol = 0.92466;
			epsInf = 1.77;
			Pvap = antoinePvap(T, 7.31549, 1794.88, -34.764);
			sigmaBulk = 4.62e-5;
			Rvdw = ScalarEOS_eval(T).vdwRadius();
			//Geometry:
			const double rOH = 1.0*Angstrom;
			const double thetaHOH = acos(-1.0/3);
			PcmSite siteO = { getAtomicNumber("O") };
			PcmSite siteH = { getAtomicNumber("H") };
			siteO.pos.push_back(vector3<>(0.,0.,0.));
			siteH.pos.push_back(vector3<>(0, -rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)));
			siteH.pos.push_back(vector3<>(0, +rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)));
			pcmSite.clear();
			pcmSite.push_back(siteO);
			pcmSite.push_back(siteH);
			break;
		}
		case CHCl3:
		{	epsBulk = 4.8069;
			Nbulk = 1.109e-3;
			pMol = 0.425;
			epsInf = 2.09;
			Pvap = antoinePvap(T, 5.96288, 1106.94, -54.598);
			sigmaBulk = 1.71e-5;
			Rvdw = TaoMasonEOS_eval(T, 536.6*Kelvin, 5328.68*KPascal, 0.216, 0.).vdwRadius();
			//Geometry:
			const double zC = 0.523*Angstrom; //distance of C from center
			const double rCCl = 1.762*Angstrom;
			const double rCH = 1.073*Angstrom;
			const double thetaHCCl = 107.98 * M_PI/180;
			PcmSite siteC = { getAtomicNumber("C") };
			PcmSite siteH = { getAtomicNumber("H") };
			PcmSite siteCl = { getAtomicNumber("Cl") };
			siteC.pos.push_back(vector3<>(0.,0.,zC));
			siteH.pos.push_back(vector3<>(0,0,zC+rCH));
			siteCl.pos.push_back(vector3<>(0, rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			siteCl.pos.push_back(vector3<>(+sqrt(0.75)*rCCl*sin(thetaHCCl), -0.5*rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			siteCl.pos.push_back(vector3<>(-sqrt(0.75)*rCCl*sin(thetaHCCl), -0.5*rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			pcmSite.clear();
			pcmSite.push_back(siteC);
			pcmSite.push_back(siteH);
			pcmSite.push_back(siteCl);
			break;
		}
		case CCl4:
		{	epsBulk = 2.238;
			Nbulk = 9.205e-4;
			pMol = 0.;
			epsInf = 2.13;
			Pvap = antoinePvap(T, 6.10445, 1265.63, -41.002);
			Rvdw = TaoMasonEOS_eval(T, 556.4*Kelvin, 4493*KPascal, 0.194, 0.).vdwRadius();
			sigmaBulk = 1.68e-5;
			//Geometry:
			const double rCCl = 1.7829*Angstrom;
			PcmSite siteC = { getAtomicNumber("C") };
			PcmSite siteCl = { getAtomicNumber("Cl") };
			siteC.pos.push_back(vector3<>(0,0,0));
			siteCl.pos.push_back(vector3<>(0,0,rCCl));
			siteCl.pos.push_back(vector3<>(0, rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			siteCl.pos.push_back(vector3<>(+sqrt(0.75)*rCCl*(sqrt(8.)/3), -0.5*rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			siteCl.pos.push_back(vector3<>(-sqrt(0.75)*rCCl*(sqrt(8.)/3), -0.5*rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			pcmSite.clear();
			pcmSite.push_back(siteC);
			pcmSite.push_back(siteCl);
			break;
		}
		case DMC:
		{	epsBulk = 3.1;
			Nbulk   = 1.059e-3;
			pMol    = 0.16;
			epsInf  = 1.87;
			Pvap = 18*mmHg;
			sigmaBulk = 2.05e-5;
			break;
		}
		case EC:
		{	epsBulk = 90.5;
			Nbulk   = 1.339e-3;
			pMol    = 2.88;
			epsInf  = 2.00;
			Pvap = antoinePvap(T, 6.05764, 1705.267, -102.261);
			sigmaBulk = 3.51e-5;
			break;
		}
		case PC:
		{	epsBulk = 64.0;
			Nbulk   = 1.039e-3;
			pMol    = 2.95;
			epsInf  = 2.02;
			Pvap = antoinePvap(T, 6.20181, 1788.900, -88.715);
			sigmaBulk = 2.88e-5;
			break;
		}
		case DMF:
		{	epsBulk = 38.0;
			Nbulk   = 1.153e-3;
			pMol    = 2.19;
			epsInf  = 2.05;
			Pvap = antoinePvap(T, 6.05286, 1400.86, -76.716);
			sigmaBulk = 2.26e-5;
			break;
		}
		case THF:
		{	epsBulk = 7.6;
			Nbulk   = 1.100e-3;
			pMol    = 0.90;
			epsInf  = 1.98;
			Pvap = antoinePvap(T, 6.12142, 1203.11, -46.795);
			sigmaBulk = 1.78e-5;
			break;
		}
	}
	
	//Set PCM fit parameters:
	switch(pcmVariant)
	{	case PCM_SGA13:
		{	nc = 7e-4;
			sigma = 0.6;
			cavityTension = 0.;
			initWarnings += "WARNING: PCM variant SGA13 is highly experimental!\n";
			break;
		}
		case PCM_GLSSA13:
		{	
			switch(solventName)
			{
				case CHCl3:
				{	nc = 2.4e-05;
					sigma = 0.6;
					cavityTension = -9.066e-6;
					break;
				}
				case CCl4:
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 1.15e-4;
							sigma = 0.6;
							cavityTension = -8.99e-06;
							break;
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with CCl4 as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
					break;
				}
				default: // For water and unparametrized fluids
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 3.7e-4;
							sigma = 0.6;
							cavityTension = 5.4e-6;
							break;
						case FluidNonlinearPCM:
							nc = 1.0e-3;
							sigma = 0.6;
							cavityTension = 9.5e-6;
							break;
						default: //Other fluids do not use these parameters
							break;
					}
					if(solventName != H2O)
					{	initWarnings +=
							"WARNING: PCM variant GLSSA has not been parametrized for this solvent; using bulk\n"
							"   surface tension as effective cavity tension and water parameters for cavity size.\n";
						cavityTension = sigmaBulk;
					}
					break;
				}
			}
			break;
		}
		case PCM_LA12:
		case PCM_PRA05:
		{	nc = 7e-4;
			sigma = 0.6;
			cavityTension = 0.;
			if(fluidType == FluidNonlinearPCM)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has not been parametrized for NonlinearPCM; using LinearPCM fit parameters.\n";
			if( (fluidType==FluidLinearPCM || fluidType==FluidNonlinearPCM) && solventName != H2O)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has been fit only for H2O; using nc and sigma from H2O fit.\n";
			break;
		}
	}
	//--- For nonlocalPCM (variant does not apply)
	if(fluidType == FluidNonlocalPCM)
	{	nc = 1.2e-3;
		sigma = sqrt(0.5);
		cavityTension = 0.;
	}
}

bool FluidSolverParams::needsVDW() const
{	switch(fluidType)
	{	case FluidNone:
			return false;
		case FluidLinearPCM:
		case FluidNonlinearPCM:
			return (pcmVariant == PCM_SGA13);
		case FluidNonlocalPCM:
		default: //All explicit fluid functionals
			return true;
	}
}

