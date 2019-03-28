/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth-Weaver

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

#include <electronic/symbols.h>
#include <fluid/FluidComponent.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_LJ.h>
#include <fluid/IdealGasMonoatomic.h>
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/IdealGasPomega.h>
#include <fluid/FluidMixture.h>

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
{	AtomicSymbol atSym;
	atomicSymbolMap.getEnum(symbol, atSym);
	return int(atSym);
}

FluidComponent::Type FluidComponent::getType(FluidComponent::Name name)
{	switch(name)
	{	case H2O:
		case CHCl3:
		case CCl4:
		case CH3CN:
		case DMC:
		case EC:
		case PC:
		case DMF:
		case THF:
		case EthylEther:
		case Chlorobenzene:
		case Isobutanol:
		case CarbonDisulfide:
		case DMSO:
		case CH2Cl2:
		case Ethanol:
		case Methanol:
		case Octanol:
		case Glyme:
		case EthyleneGlycol:
		case CustomSolvent:
			return Solvent;
		case Sodium:
		case HydratedSodium:
		case Potassium:
		case HydratedPotassium:
		case Hydronium:
		case HydratedHydronium: 
		case CustomCation:
			return Cation;
		case Chloride:
		case Fluoride: 
		case Perchlorate: 
		case Hydroxide:
		case HydratedHydroxide:
		case CustomAnion: 
			return Anion;
		default:
			assert(!"Unknown component type");
			return Solvent;
	}
}

std::vector<complex> FluidComponent::getChiPrefactor(const std::vector<complex>& omegaArr, double chi0nuc, double chi0el) const
{	std::vector<complex> result; result.reserve(omegaArr.size());
	for(const complex& omega: omegaArr)
	{	complex chi = complex(chi0nuc)/(1. - complex(0,tauNuc)*omega); //nuclear part
		//Add electronic poles
		for(const PoleLD& pole: polesEl)
			chi += complex(chi0el * pole.A0 * pole.omega0*pole.omega0) /
				(pole.omega0*pole.omega0 - omega*(omega + complex(0,pole.gamma0)));
		result.push_back(chi);
	}
	return result;
}

double FluidComponent::pureNbulk(double T) const
{	if(type == Solvent)
	{	switch(name) //TODO: add temperature dependence
		{	case H2O: return 4.9383e-3;
			case CHCl3: return 1.109e-3;
			case CCl4: return 9.205e-4;
			case CH3CN: return 1.709e-3;
			case DMC: return 1.059e-3;
			case EC: return 1.339e-3;
			case PC: return 1.039e-3;
			case DMF: return 1.153e-3;
			case THF: return 1.100e-3;
			case EthylEther: return 8.5e-4;
			case Isobutanol: return 9.668e-4;
			case Chlorobenzene: return 8.74e-4;
			case CarbonDisulfide: return 1.48e-3;
			case DMSO: return 1.256e-3;
			case CH2Cl2: return 1.392e-3;
			case Ethanol: return 1.528e-3;
			case Methanol: return 2.203e-3;
			case Octanol: return 5.646e-4;
			case Glyme: return 8.586e-4;
			case EthyleneGlycol: return 1.60e-3;
			default: throw string("Not yet implemented.");
		}
	}
	else return 1.*mol/liter; //ions
}


FluidComponent::FluidComponent(FluidComponent::Name name, double T, FluidComponent::Functional functional)
: name(name), type(getType(name)), functional(functional), epsLJ(0.), representation(MuEps),
s2quadType(Quad7design_24), quad_nBeta(0), quad_nAlpha(0), quad_nGamma(0), translationMode(LinearSpline),
epsBulk(1.), Nbulk(pureNbulk(T)), pMol(0.), epsInf(1.), Pvap(0.), sigmaBulk(0.), Rvdw(0.), Res(0.),
tauNuc(8.3e+3*fs), Nnorm(0), quad(0), trans(0), idealGas(0), fex(0), offsetIndep(0), offsetDensity(0)
{
	//Nuclear widths = (1./6) vdW radius
	const double sigmaNucH = (1./6) * 1.20*Angstrom;
	const double sigmaNucC = (1./6) * 1.70*Angstrom;
	const double sigmaNucN = (1./6) * 1.55*Angstrom;
	const double sigmaNucO = (1./6) * 1.52*Angstrom;
	const double sigmaNucCl = (1./6) * 1.75*Angstrom;

	//Set physical parameters (in atomic units) decsribing solvent:
	switch(name)
	{	case H2O:
		{	epsBulk = 78.4;
			pMol = 0.92466;
			epsInf = 1.77;
			Pvap = antoinePvap(T, 7.31549, 1794.88, -34.764);
			sigmaBulk = 4.62e-5;
			eos = std::make_shared<JeffereyAustinEOS>(T, 2*(1.36*Angstrom));
			Rvdw = 1.385*Angstrom;
			Res = 1.42;
			//Simgle-pole frequency model:
			PoleLD pole;
			pole.omega0 = 15.*eV;
			pole.gamma0 = 7.*eV;
			pole.A0 = 1.;
			polesEl.assign(1, pole);
			//Site properties:
			molecule.name = "H2O";
			auto siteO = std::make_shared<Molecule::Site>("O",int(AtomicSymbol::O));
				siteO->Znuc = 6.; siteO->sigmaNuc = sigmaNucO;
				siteO->Zelec = 6.826; //siteO->aElec = 0.32;
				siteO->aElec = 0.36935; siteO->sigmaElec = 0.36823; siteO->rcElec = 0.51523;
				siteO->alpha = 3.73; siteO->aPol = 0.32;
				siteO->Zsite = siteO->Zelec-siteO->Znuc;
			molecule.sites.push_back(siteO);
			auto siteH = std::make_shared<Molecule::Site>("H",int(AtomicSymbol::H));
				siteH->Znuc = 1.; siteH->sigmaNuc = sigmaNucH;
				siteH->Zelec = 0.587; //siteH->aElec = 0.31;
				siteH->aElec = 0.34641; siteH->sigmaElec = 2.*(siteH->aElec)/sqrt(M_PI); siteH->rcElec = 0.0;
				siteH->alpha = 3.30; siteH->aPol = 0.39;
				siteH->Zsite = siteH->Zelec-siteH->Znuc;
			molecule.sites.push_back(siteH);
			//Geometry:
			const double rOH = 0.967*Angstrom;
			const double thetaHOH = 104.2 * M_PI/180;
			siteO->positions.push_back(vector3<>(0.,0.,0.));
			siteH->positions.push_back(vector3<>(0, -rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)));
			siteH->positions.push_back(vector3<>(0, +rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)));
			//Functional dependent options:
			switch(functional)
			{	case FittedCorrelations:
					break;
				case ScalarEOS:
					siteO->Rhs = 0.5*eos->sigmaEOS;
					break;
				case BondedVoids: 
				{	siteO->Rhs = Fex_H2O_BondedVoids::RO;
					//Add void sites:
					auto siteV = std::make_shared<Molecule::Site>("V");
					molecule.sites.push_back(siteV);
					siteV->Rhs = Fex_H2O_BondedVoids::RV0 * exp(-T/Fex_H2O_BondedVoids::TV);
					const double rOV = siteO->Rhs + siteV->Rhs;
					const double thetaVOV = acos(-1./3);
					siteV->positions.push_back(vector3<>(-rOV*sin(0.5*thetaVOV), 0, -rOV*cos(0.5*thetaVOV)));
					siteV->positions.push_back(vector3<>(+rOV*sin(0.5*thetaVOV), 0, -rOV*cos(0.5*thetaVOV)));
					break;
				}
				default:
					die("Unsupported excess functional for water.\n");
			}
			break;
		}
		case CHCl3:
		{	epsBulk = 4.8069;
			pMol = 0.49091;
			epsInf = 2.09;
			Pvap = antoinePvap(T, 5.96288, 1106.94, -54.598);
			sigmaBulk = 1.71e-5;
			eos = std::make_shared<TaoMasonEOS>(T, 536.6*Kelvin, 5328.68*KPascal, 0.216, 2.70*Angstrom);
			Rvdw = 2.53*Angstrom;
			Res = 2.22;
			//Site properties:
			molecule.name = "CHCl3";
			auto siteCenter = std::make_shared<Molecule::Site>("center",0);
				siteCenter->Rhs = Rvdw;
			molecule.sites.push_back(siteCenter);
			auto siteC = std::make_shared<Molecule::Site>("C",int(AtomicSymbol::C));
				siteC->Znuc = 4.; siteC->sigmaNuc = sigmaNucC;
				siteC->Zelec = 4.256; //siteC->aElec = 0.43;
				siteC->aElec = 0.49; siteC->sigmaElec = 0.48; siteC->rcElec = 0.67;
				siteC->alpha = 6.05; siteC->aPol = 0.36;
				siteC->Zsite = siteC->Zelec-siteC->Znuc;
			molecule.sites.push_back(siteC);
			auto siteH = std::make_shared<Molecule::Site>("H",int(AtomicSymbol::H));
				siteH->Znuc = 1.; siteH->sigmaNuc = sigmaNucH;
				siteH->Zelec = 0.756; //siteH->aElec = 0.26;
				siteH->aElec = 0.36; siteH->sigmaElec = 2.*(siteH->aElec)/sqrt(M_PI); siteH->rcElec = 0.0;
				siteH->alpha = 9.13; siteH->aPol = 0.41;
				siteH->Zsite = siteH->Zelec-siteH->Znuc;
			molecule.sites.push_back(siteH);
			auto siteCl = std::make_shared<Molecule::Site>("Cl",int(AtomicSymbol::Cl));
				siteCl->Znuc = 7.; siteCl->sigmaNuc = sigmaNucCl;
				siteCl->Zelec = 6.996; //siteCl->aElec = 0.44;
				siteCl->aElec = 0.45; siteCl->sigmaElec = 0.51; siteCl->rcElec = 1.01;
				siteCl->alpha = 15.8; siteCl->aPol = 0.46;
				siteCl->Zsite = siteCl->Zelec-siteCl->Znuc;
			molecule.sites.push_back(siteCl);
			//Geometry:
			const double zC = 0.523*Angstrom; //distance of C from center
			const double rCCl = 1.804*Angstrom;
			const double rCH = 1.091*Angstrom;
			const double thetaHCCl = 107.8 * M_PI/180;
			siteCenter->positions.push_back(vector3<>(0.,0.,0.));
			siteC->positions.push_back(vector3<>(0.,0.,zC));
			siteH->positions.push_back(vector3<>(0,0,zC+rCH));
			siteCl->positions.push_back(vector3<>(0, rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			siteCl->positions.push_back(vector3<>(+sqrt(0.75)*rCCl*sin(thetaHCCl), -0.5*rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			siteCl->positions.push_back(vector3<>(-sqrt(0.75)*rCCl*sin(thetaHCCl), -0.5*rCCl*sin(thetaHCCl), zC+rCCl*cos(thetaHCCl)));
			break;
		}
		case CCl4:
		{	epsBulk = 2.238;
			pMol = 0.;
			epsInf = 2.13;
			Pvap = antoinePvap(T, 6.10445, 1265.63, -41.002);
			sigmaBulk = 1.68e-5;
			eos = std::make_shared<TaoMasonEOS>(T, 556.4*Kelvin, 4493*KPascal, 0.194, 2.78*Angstrom);
			Rvdw = 2.69*Angstrom;
			Res = 1.90;
			//Site properties:
			molecule.name = "CCl4";
			auto siteC = std::make_shared<Molecule::Site>("C",int(AtomicSymbol::C));
				siteC->Znuc = 4.; siteC->sigmaNuc = sigmaNucC;
				siteC->Zelec = 4.980; //siteC->aElec = 0.44;
				siteC->aElec = 0.61; siteC->sigmaElec = 0.37; siteC->rcElec = 0.53;
				siteC->alpha = 5.24; siteC->aPol = 0.35;
				siteC->Rhs = Rvdw;
				siteC->Zsite = siteC->Zelec-siteC->Znuc;
			molecule.sites.push_back(siteC);
			auto siteCl = std::make_shared<Molecule::Site>("Cl",int(AtomicSymbol::Cl));
				siteCl->Znuc = 7.; siteCl->sigmaNuc = sigmaNucCl;
				siteCl->Zelec = 6.755; //siteCl->aElec = 0.44;
				siteCl->aElec = 0.44; siteCl->sigmaElec = 0.52; siteCl->rcElec = 1.04;
				siteCl->alpha = 18.1; siteCl->aPol = 0.47;
				siteCl->Zsite = siteCl->Zelec-siteCl->Znuc;
			molecule.sites.push_back(siteCl);
			//Geometry:
			const double rCCl = 1.801*Angstrom;
			siteC->positions.push_back(vector3<>(0,0,0));
			siteCl->positions.push_back(vector3<>(0,0,rCCl));
			siteCl->positions.push_back(vector3<>(0, rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			siteCl->positions.push_back(vector3<>(+sqrt(0.75)*rCCl*(sqrt(8.)/3), -0.5*rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			siteCl->positions.push_back(vector3<>(-sqrt(0.75)*rCCl*(sqrt(8.)/3), -0.5*rCCl*(sqrt(8.)/3), rCCl*(-1./3)));
			break;
		}
		case CH3CN:
		{	epsBulk = 38.8;
			pMol = 1.89;
			epsInf = 1.81;
			Pvap = antoinePvap(T, 6.52111, 1492.375, -24.208); //logPrintf("selfSol = %lg\n", T*log(Pvap/(Nbulk*T)));
			sigmaBulk = 1.88e-5;
			eos = std::make_shared<TaoMasonEOS>(T, 545.5*Kelvin, 4830*KPascal, 0.278, 3.5*Angstrom); //TODO: Tao-Mason EOS not appropriate for this fluid
			Rvdw = 2.12*Angstrom;
			Res = 2.6; //empirical value, SaLSA predicts 3.13
			//Site properties:
			molecule.name = "CH3CN";
			auto siteCenter = std::make_shared<Molecule::Site>("center",0);
				siteCenter->Rhs = 1.12*Angstrom;
			molecule.sites.push_back(siteCenter);
			auto siteC1 = std::make_shared<Molecule::Site>("C1",int(AtomicSymbol::C)); //methyl carbon
				siteC1->Znuc = 4.; siteC1->sigmaNuc = sigmaNucC;
				siteC1->Zelec = 4.7128; siteC1->aElec = 0.44;
				siteC1->alpha = 4.49; siteC1->aPol = 0.35;
				siteC1->Zsite = siteC1->Zelec-siteC1->Znuc;
			molecule.sites.push_back(siteC1);
			auto siteC2 = std::make_shared<Molecule::Site>("C2",int(AtomicSymbol::C)); //nitrile carbon
				siteC2->Znuc = 4.; siteC2->sigmaNuc = sigmaNucC;
				siteC2->Zelec = 3.4832; siteC2->aElec = 0.39;
				siteC2->alpha = 7.18; siteC2->aPol = 0.39;
				siteC2->Zsite = siteC2->Zelec-siteC2->Znuc;
			molecule.sites.push_back(siteC2);
			auto siteH = std::make_shared<Molecule::Site>("H",int(AtomicSymbol::H));
				siteH->Znuc = 1.; siteH->sigmaNuc = sigmaNucH;
				siteH->Zelec = 0.7659; siteH->aElec = 0.28;
				siteH->alpha = 4.33; siteH->aPol = 0.37;
				siteH->Zsite = siteH->Zelec-siteH->Znuc;
			molecule.sites.push_back(siteH);
			auto siteN = std::make_shared<Molecule::Site>("N",int(AtomicSymbol::N));
				siteN->Znuc = 5.; siteN->sigmaNuc = sigmaNucN;
				siteN->Zelec = 5.5063; siteN->aElec = 0.37;
				siteN->alpha = 5.85; siteN->aPol = 0.35;
				siteN->Zsite = siteN->Zelec-siteN->Znuc;
			molecule.sites.push_back(siteN);
			//Geometry:
			const double zC2 = 0.165*Angstrom; //distance of nitrile carbon from center
			const double rCC = 1.462*Angstrom;
			const double rCN = 1.161*Angstrom;
			const double rCH = 1.098*Angstrom;
			const double thetaCCH = 110.22*M_PI/180;
			siteCenter->positions.push_back(vector3<>(0.,0.,0.));
			siteC2->positions.push_back(vector3<>(0.,0.,zC2));
			siteC1->positions.push_back(vector3<>(0.,0.,zC2-rCC));
			siteN->positions.push_back(vector3<>(0,0,zC2+rCN));
			siteH->positions.push_back(vector3<>(0, rCH*sin(thetaCCH), zC2-rCC+rCH*cos(thetaCCH)));
			siteH->positions.push_back(vector3<>(+sqrt(0.75)*rCH*sin(thetaCCH), -0.5*rCH*sin(thetaCCH), zC2-rCC+rCH*cos(thetaCCH)));
			siteH->positions.push_back(vector3<>(-sqrt(0.75)*rCH*sin(thetaCCH), -0.5*rCH*sin(thetaCCH), zC2-rCC+rCH*cos(thetaCCH)));
			break;
		}
		case CH2Cl2:
		{	epsBulk = 9.08;
			pMol    = 0.89;
			epsInf  = 1.424;
			sigmaBulk = 1.70e-5;
			molecule.name = "CH2Cl2";
			break;
		}
		case Ethanol:
		{	epsBulk = 24.3;
			pMol    = 0.76;
			epsInf  = 1.361;
			sigmaBulk = 1.44e-5;
			molecule.name = "Ethanol";
			break;
		}
		case Methanol:
		{	epsBulk = 32.66;
			pMol = 0.791;
			epsInf = 1.328;
			sigmaBulk = 1.445e-5;
			molecule.name = "Methanol";
			break;
		}
		case Octanol:
		{	epsBulk = 10.30;
			pMol = 0.661;
			epsInf = 2.036;
			sigmaBulk = 1.766e-5;
			Pvap = antoinePvap(T, 8.47682, 2603.359, -48.799);
			Rvdw = 3.348*Angstrom;
			molecule.name = "Octanol";
			break;
		}
		case DMC:
		{	epsBulk = 3.1;
			pMol    = 0.16;
			epsInf  = 1.87;
			Pvap = 18*mmHg;
			sigmaBulk = 2.05e-5;
			molecule.name = "DMC";
			break;
		}
		case EC:
		{	epsBulk = 90.5;
			pMol    = 2.88;
			epsInf  = 2.00;
			Pvap = antoinePvap(T, 6.05764, 1705.267, -102.261);
			sigmaBulk = 3.51e-5;
			molecule.name = "EC";
			break;
		}
		case PC:
		{	epsBulk = 64.0;
			pMol    = 2.95;
			epsInf  = 2.02;
			Pvap = antoinePvap(T, 6.20181, 1788.900, -88.715);
			sigmaBulk = 2.88e-5;
			molecule.name = "PC";
			break;
		}
		case DMF:
		{	epsBulk = 38.0;
			pMol    = 2.19;
			epsInf  = 2.05;
			Pvap = antoinePvap(T, 6.05286, 1400.86, -76.716);
			sigmaBulk = 2.26e-5;
			molecule.name = "DMF";
			break;
		}
		case THF:
		{	epsBulk = 7.6;
			pMol    = 0.90;
			epsInf  = 1.98;
			Pvap = antoinePvap(T, 6.12142, 1203.11, -46.795);
			sigmaBulk = 1.78e-5;
			molecule.name = "THF";
			break;
		}
		case DMSO:
		{	epsBulk = 48;
			pMol    = 1.56;
			epsInf  = 2.19;
			sigmaBulk = 2.80e-5;
			Pvap = antoinePvap(T, 7.23039, 2239.161, -29.215);
			Rvdw = 2.378*Angstrom;
			molecule.name = "DMSO";
			break;
		}
		case EthylEther:
		{	epsBulk = 4.34;
			pMol = 0.487;
			epsInf = 1.82;
			Pvap = antoinePvap(T, 6.96559, 1071.54, 227.774);
			sigmaBulk = 1.092e-5;
			molecule.name = "EthylEther";
			break;
		}
		case Chlorobenzene:
		{	epsBulk = 5.69;
			pMol = 0.72;
			epsInf = 2.32;
			Pvap = antoinePvap(T, 4.11083, 1435.675, -55.124);
			sigmaBulk = 2.1e-5;
			molecule.name = "Chlorobenzene";
			break;
		}
		case Isobutanol:
		{	epsBulk = 17.93;
			pMol = 0.646;
			epsInf = 1.949;
			sigmaBulk = 1.445e-5;
			molecule.name = "Isobutanol";
			break;
		}
		case CarbonDisulfide:
		{	epsBulk = 2.641;
			epsInf = 2.641;
			pMol = 0.;
			molecule.name = "CarbonDisulfide";
			break;
		}
		case Glyme:
		{	epsBulk = 7.20;
			epsInf = 1.90;
			pMol = 0.;
			molecule.name = "Glyme";
			break;
		}
		case EthyleneGlycol:
		{	epsBulk = 41.4;
			epsInf = 1.43;
			pMol = 0.;
			molecule.name = "EthyleneGlycol";
			break;
		}
		case Sodium: 
		{	Rvdw = 1.16*Angstrom;
			//Site properties:
			molecule.name = "Na+";
			auto siteNa = std::make_shared<Molecule::Site>("Na",int(AtomicSymbol::Na));
			siteNa->Znuc = 9.1383; siteNa->sigmaNuc = (1./6)*Rvdw;
			siteNa->Zelec = 8.1383; //siteNa->aElec = 0.206;
			siteNa->aElec = 0.19682; siteNa->sigmaElec = 0.41314; siteNa->rcElec = 0.71491;
			siteNa->Zsite = siteNa->Zelec-siteNa->Znuc;
			siteNa->Rhs = 0.986*Angstrom; //reproduces location of first peak in Na-O radial distribution functions
			molecule.sites.push_back(siteNa);
			//Geometry:
			siteNa->positions.push_back(vector3<>(0,0,0));
			break;
		}
		case Potassium:
		{	Rvdw = 1.51*Angstrom;
			//Site properties:
			molecule.name = "K+";
			auto siteK = std::make_shared<Molecule::Site>("K",int(AtomicSymbol::K));
			siteK->Znuc = 9.3230; siteK->sigmaNuc = (1./6)*Rvdw;
			siteK->Zelec = 8.3230; siteK->aElec = 0.35529;
			siteK->sigmaElec = 0.45583; siteK->rcElec =  0.82179;
			siteK->Zsite = siteK->Zelec-siteK->Znuc;
			//siteK->Rhs = Rvdw; //obvious default
			siteK->Rhs = 1.4330*Angstrom; //reproduces location of first peak in K-O radial distribution functions (Chem. Rev. 88, 1988)
			molecule.sites.push_back(siteK);
			//Geometry:
			siteK->positions.push_back(vector3<>(0,0,0));
			break;
		}
		case Chloride:
		{	Rvdw = 1.67*Angstrom;
			//Site properties:
			molecule.name = "Cl-";
			auto siteCl = std::make_shared<Molecule::Site>("Cl",int(AtomicSymbol::Cl));
			siteCl->Znuc = 8.6677; siteCl->sigmaNuc = (1./6)*Rvdw;
			siteCl->Zelec = 9.6677; //siteCl->aElec = 0.438;
			siteCl->aElec = 0.54615; siteCl->sigmaElec =  2.*(siteCl->aElec)/sqrt(M_PI); siteCl->rcElec = 0.0;			
			siteCl->Zsite = siteCl->Zelec-siteCl->Znuc;
			siteCl->Rhs = 1.837*Angstrom; //Reproduces first peak in the Cl-O radial distribution functions
			molecule.sites.push_back(siteCl);
			//Geometry:
			siteCl->positions.push_back(vector3<>(0,0,0));
			break;
		}
		case Fluoride:
		{	Rvdw = 1.19*Angstrom;
			//Site properties:
			molecule.name = "F-";
			auto siteF = std::make_shared<Molecule::Site>("F",int(AtomicSymbol::F));
			siteF->Znuc = 7.; siteF->sigmaNuc = (1./6)*Rvdw;
			siteF->Zelec = 8.; siteF->aElec = 0.38886;
			siteF->sigmaElec =  2.*(siteF->aElec)/sqrt(M_PI); siteF->rcElec = 0.0;		  
			siteF->Zsite = siteF->Zelec-siteF->Znuc;
			siteF->Rhs = 1.27*Angstrom; //Reproduces first peak in the F-O radial distribution functions
			molecule.sites.push_back(siteF);
			//Geometry:
			siteF->positions.push_back(vector3<>(0,0,0));
			break;
		}
		case Perchlorate:
		{	Rvdw = 2.41*Angstrom; //very close to Pauling Ionic Radius
			//Site properties:
			molecule.name = "ClO4-";
			auto siteCl = std::make_shared<Molecule::Site>("Cl",int(AtomicSymbol::Cl));
			siteCl->Znuc = 7.; siteCl->sigmaNuc = sigmaNucCl;
			siteCl->Zelec =  6.2308; siteCl->aElec = 0.48793;
			siteCl->sigmaElec = 0.45657; siteCl->rcElec = 0.93371;
			siteCl->Zsite = siteCl->Zelec-siteCl->Znuc;   
			siteCl->Rhs = 2.34*Angstrom;
			molecule.sites.push_back(siteCl);
			auto siteO = std::make_shared<Molecule::Site>("O",int(AtomicSymbol::O));
			siteO->Znuc = 6.; siteO->sigmaNuc = sigmaNucO;
			siteO->Zelec =  6.4423; siteO->aElec = 0.36013;
			siteO->sigmaElec = 0.37371; siteO->rcElec = 0.51636;
			siteO->Zsite = siteO->Zelec-siteO->Znuc;
			molecule.sites.push_back(siteO);
			//geometry
			double rClO = 1.563*Angstrom;
			siteCl->positions.push_back(vector3<>(0,0,0));
			siteO->positions.push_back(vector3<>(rClO,rClO,rClO)/sqrt(3));
			siteO->positions.push_back(vector3<>(rClO,-rClO,-rClO)/sqrt(3));
			siteO->positions.push_back(vector3<>(-rClO,rClO,-rClO)/sqrt(3));
			siteO->positions.push_back(vector3<>(-rClO,-rClO,rClO)/sqrt(3));
			break;
		}
		default:
			throw string("Not yet implemented.");
	}
}

void FluidComponent::addToFluidMixture(FluidMixture* fluidMixture)
{	assert(!idealGas);	
        const GridInfo& gInfo = fluidMixture->gInfo;
	if(!molecule) molecule.setup(gInfo, Rvdw);
	
	//Setup ideal gas:
	if(molecule.isMonoatomic())
	{	idealGas = std::make_shared<IdealGasMonoatomic>(fluidMixture, this);
	}
	else
	{	quad = std::make_shared<SO3quad>(s2quadType, molecule, quad_nBeta, quad_nAlpha, quad_nGamma);
		switch(translationMode)
		{	case LinearSpline: trans = std::make_shared<TranslationOperatorSpline>(gInfo, TranslationOperatorSpline::Linear); break;
			case ConstantSpline: trans = std::make_shared<TranslationOperatorSpline>(gInfo, TranslationOperatorSpline::Constant); break;
			case Fourier: trans = std::make_shared<TranslationOperatorFourier>(gInfo); break;
		}
		switch(representation)
		{	case PsiAlpha: idealGas = std::make_shared<IdealGasPsiAlpha>(fluidMixture, this, *quad, *trans); break;
			case Pomega: idealGas = std::make_shared<IdealGasPomega>(fluidMixture, this, *quad, *trans); break;
			case MuEps: idealGas = std::make_shared<IdealGasMuEps>(fluidMixture, this, *quad, *trans); break;
		}
	}
	
	//Initialize excess functional:
	switch(functional)
	{	case ScalarEOS:
			assert(eos);
			fex = std::make_shared<Fex_ScalarEOS>(fluidMixture, this, *eos);
			break;
		case BondedVoids:
			assert(name == H2O);
			fex = std::make_shared<Fex_H2O_BondedVoids>(fluidMixture, this);
			break;
		case FittedCorrelations:
			assert(name == H2O);
			fex = std::make_shared<Fex_H2O_FittedCorrelations>(fluidMixture, this);
			break;
		case MeanFieldLJ:
			assert(molecule.sites[0]->Rhs > 0.);
			fex = std::make_shared<Fex_LJ>(fluidMixture, this, epsLJ);
			break;
		case FunctionalNone:
			//No excess functional, or manually created functional (not managed by FluidComponent)
			break;
	}

	fluidMixture->addComponent(this);
}
