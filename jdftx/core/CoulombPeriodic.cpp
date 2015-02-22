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

#include <core/CoulombPeriodic.h>
#include <core/Coulomb_internal.h>
#include <core/CoulombKernel.h>
#include <core/BlasExtra.h>

//! Standard 3D Ewald sum
class EwaldPeriodic : public Ewald
{
	matrix3<> R, G, RTR, GGT; //!< Lattice vectors, reciprocal lattice vectors and corresponding metrics
	double sigma; //!< gaussian width for Ewald sums
	vector3<int> Nreal; //!< max unit cell indices for real-space sum
	vector3<int> Nrecip; //!< max unit cell indices for reciprocal-space sum

public:
	EwaldPeriodic(const matrix3<>& R, int nAtoms)
	: R(R), G((2*M_PI)*inv(R)), RTR((~R)*R), GGT(G*(~G))
	{	logPrintf("\n---------- Setting up ewald sum ----------\n");
		//Determine optimum gaussian width for Ewald sums:
		// From below, the number of reciprocal cells ~ Prod_k |R.column[k]|
		//    and number of real space cells ~ Prod_k |G.row[k]|
		// including the fact that the real space cost ~ Natoms^2/cell
		//    and the reciprocal space cost ~ Natoms/cell
		sigma = 1.;
		for(int k=0; k<3; k++)
			sigma *= R.column(k).length() / G.row(k).length();
		sigma = pow(sigma/std::max(1,nAtoms), 1./6);
		logPrintf("Optimum gaussian width for ewald sums = %lf bohr.\n", sigma);
		
		//Carry real space sums to Rmax = 10 sigma and Gmax = 10/sigma
		//This leads to relative errors ~ 1e-22 in both sums, well within double precision limits
		for(int k=0; k<3; k++)
		{	Nreal[k] = 1+ceil(CoulombKernel::nSigmasPerWidth * G.row(k).length() * sigma / (2*M_PI));
			Nrecip[k] = 1+ceil(CoulombKernel::nSigmasPerWidth * R.column(k).length() / (2*M_PI*sigma));
		}
		logPrintf("Real space sum over %d unit cells with max indices ", (2*Nreal[0]+1)*(2*Nreal[1]+1)*(2*Nreal[2]+1));
		Nreal.print(globalLog, " %d ");
		logPrintf("Reciprocal space sum over %d terms with max indices ", (2*Nrecip[0]+1)*(2*Nrecip[1]+1)*(2*Nrecip[2]+1));
		Nrecip.print(globalLog, " %d ");
	}

	double energyAndGrad(std::vector<Atom>& atoms) const
	{	double eta = sqrt(0.5)/sigma, etaSq=eta*eta;
		double sigmaSq = sigma * sigma;
		double detR = fabs(det(R)); //cell volume
		//Position independent terms:
		double Ztot = 0., ZsqTot = 0.;
		for(const Atom& a: atoms)
		{	Ztot += a.Z;
			ZsqTot += a.Z * a.Z;
		}
		double E
			= 0.5 * 4*M_PI * Ztot*Ztot * (-0.5*sigmaSq) / detR //G=0 correction
			- 0.5 * ZsqTot * eta * (2./sqrt(M_PI)); //Self-energy correction
		//Reduce positions to first centered unit cell:
		for(Atom& a: atoms)
			for(int k=0; k<3; k++)
				a.pos[k] -= floor(0.5 + a.pos[k]);
		//Real space sum:
		vector3<int> iR; //integer cell number
		for(const Atom& a2: atoms)
			for(Atom& a1: atoms)
				for(iR[0]=-Nreal[0]; iR[0]<=Nreal[0]; iR[0]++)
					for(iR[1]=-Nreal[1]; iR[1]<=Nreal[1]; iR[1]++)
						for(iR[2]=-Nreal[2]; iR[2]<=Nreal[2]; iR[2]++)
						{	vector3<> x = iR + (a1.pos - a2.pos);
							double rSq = RTR.metric_length_squared(x);
							if(!rSq) continue; //exclude self-interaction
							double r = sqrt(rSq);
							E += 0.5 * a1.Z * a2.Z * erfc(eta*r)/r;
							a1.force += (RTR * x) *
								(a1.Z * a2.Z * (erfc(eta*r)/r + (2./sqrt(M_PI))*eta*exp(-etaSq*rSq))/rSq);
						}
		//Reciprocal space sum:
		vector3<int> iG; //integer reciprocal cell number
		for(iG[0]=-Nrecip[0]; iG[0]<=Nrecip[0]; iG[0]++)
			for(iG[1]=-Nrecip[1]; iG[1]<=Nrecip[1]; iG[1]++)
				for(iG[2]=-Nrecip[2]; iG[2]<=Nrecip[2]; iG[2]++)
				{	double Gsq = GGT.metric_length_squared(iG);
					if(!Gsq) continue; //skip G=0
					//Compute structure factor:
					complex SG = 0.;
					for(const Atom& a: atoms)
						SG += a.Z * cis(-2*M_PI*dot(iG,a.pos));
					//Accumulate energy:
					double eG = 4*M_PI * exp(-0.5*sigmaSq*Gsq)/(Gsq * detR);
					E += 0.5 * eG * SG.norm();
					//Accumulate forces:
					for(Atom& a: atoms)
						a.force -= (eG * a.Z * 2*M_PI * (SG.conj() * cis(-2*M_PI*dot(iG,a.pos))).imag()) * iG;
				}
		return E;
	}
};


//------------- class CoulombPeriodic ---------------

CoulombPeriodic::CoulombPeriodic(const GridInfo& gInfoOrig, const CoulombParams& params)
: Coulomb(gInfoOrig, params)
{
	initExchangeEval();
}

ScalarFieldTilde CoulombPeriodic::apply(ScalarFieldTilde&& in) const
{	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombPeriodic_calc(), in->dataPref(false));
	return in;
}

std::shared_ptr<Ewald> CoulombPeriodic::createEwald(matrix3<> R, size_t nAtoms) const
{	return std::make_shared<EwaldPeriodic>(R, nAtoms);
}
