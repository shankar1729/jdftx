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
#include <core/BlasExtra.h>

//! Standard 3D Ewald sum
struct EwaldPeriodic
{
	const GridInfo& gInfo;
	double sigma; //!< gaussian width for Ewald sums
	vector3<int> Nreal; //!< max unit cell indices for real-space sum
	vector3<int> Nrecip; //!< max unit cell indices for reciprocal-space sum
	
	EwaldPeriodic(const GridInfo& gInfo, int nAtoms) : gInfo(gInfo)
	{	logPrintf("\n---------- Setting up ewald sum ----------\n");
		//Determine optimum gaussian width for Ewald sums:
		// From below, the number of reciprocal cells ~ Prod_k |R.column[k]|
		//    and number of real space cells ~ Prod_k |G.row[k]|
		// including the fact that the real space cost ~ Natoms^2/cell
		//    and the reciprocal space cost ~ Natoms/cell
		sigma = 1.;
		for(int k=0; k<3; k++)
			sigma *= gInfo.R.column(k).length() / gInfo.G.row(k).length();
		sigma = pow(sigma/std::max(1,nAtoms), 1./6);
		logPrintf("Optimum gaussian width for ewald sums = %lf bohr.\n", sigma);
		
		//Carry real space sums to Rmax = 10 sigma and Gmax = 10/sigma
		//This leads to relative errors ~ 1e-22 in both sums, well within double precision limits
		for(int k=0; k<3; k++)
		{	Nreal[k] = 1+ceil(10. * gInfo.G.row(k).length() * sigma / (2*M_PI));
			Nrecip[k] = 1+ceil(10. * gInfo.R.column(k).length() / (2*M_PI*sigma));
		}
		logPrintf("Real space sums over %d unit cells with max indices ", 8*Nreal[0]*Nreal[1]*Nreal[2]);
		Nreal.print(globalLog, " %d ");
		logPrintf("Reciprocal space sums over %d unit cells with max indices ", 8*Nrecip[0]*Nrecip[1]*Nrecip[2]);
		Nrecip.print(globalLog, " %d ");
	}
	
	double energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
	{	double eta = sqrt(0.5)/sigma, etaSq=eta*eta;
		double sigmaSq = sigma * sigma;
		//Position independent terms:
		double Ztot = 0., ZsqTot = 0.;
		for(const Coulomb::PointCharge& pc: pointCharges)
		{	Ztot += pc.Z;
			ZsqTot += pc.Z * pc.Z;
		}
		double E
			= 0.5 * 4*M_PI * Ztot*Ztot * (-0.5*sigmaSq) / gInfo.detR //G=0 correction
			- 0.5 * ZsqTot * eta * (2./sqrt(M_PI)); //Self-energy correction
		//Reduce positions to first unit cell:
		for(Coulomb::PointCharge& pc: pointCharges)
			for(int k=0; k<3; k++)
				pc.pos[k] -= floor(pc.pos[k]);
		//Real space sum:
		vector3<int> iR; //integer cell number
		for(iR[0]=-Nreal[0]+1; iR[0]<=Nreal[0]; iR[0]++)
			for(iR[1]=-Nreal[1]+1; iR[1]<=Nreal[1]; iR[1]++)
				for(iR[2]=-Nreal[2]+1; iR[2]<=Nreal[2]; iR[2]++)
					for(const Coulomb::PointCharge& pc2: pointCharges)
						for(Coulomb::PointCharge& pc1: pointCharges)
						{	vector3<> x = iR + pc1.pos - pc2.pos;
							double rSq = gInfo.RTR.metric_length_squared(x);
							double r = sqrt(rSq);
							if(!r) continue; //exclude self-interaction
							E += 0.5 * pc1.Z * pc2.Z * erfc(eta*r)/r;
							pc1.force += (gInfo.RTR * x) *
								(pc1.Z * pc2.Z * (erfc(eta*r)/r + (2./sqrt(M_PI))*eta*exp(-etaSq*rSq))/rSq);
						}
		//Reciprocal space sum:
		vector3<int> iG; //integer reciprocal cell number
		for(iG[0]=-Nrecip[0]+1; iG[0]<=Nrecip[0]; iG[0]++)
			for(iG[1]=-Nrecip[1]+1; iG[1]<=Nrecip[1]; iG[1]++)
				for(iG[2]=-Nrecip[2]+1; iG[2]<=Nrecip[2]; iG[2]++)
				{	double Gsq = gInfo.GGT.metric_length_squared(iG);
					if(!Gsq) continue; //skip G=0
					//Compute structure factor:
					complex SG = 0.;
					for(const Coulomb::PointCharge& pc: pointCharges)
						SG += pc.Z * cis(-2*M_PI*dot(iG,pc.pos));
					//Accumulate energy:
					double eG = 4*M_PI * exp(-0.5*sigmaSq*Gsq)/(Gsq * gInfo.detR);
					E += 0.5 * eG * SG.norm();
					//Accumulate forces:
					for(Coulomb::PointCharge& pc: pointCharges)
						pc.force -= (eG * pc.Z * 2*M_PI * (SG.conj() * cis(-2*M_PI*dot(iG,pc.pos))).imag()) * iG;
				}
		return E;
	}
};

CoulombPeriodic::CoulombPeriodic(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params)
{
}

DataGptr CoulombPeriodic::operator()(DataGptr&& in) const
{	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombPeriodic_calc(), in->dataPref(false));
	return in;
}

double CoulombPeriodic::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	if(!ewald)
		((CoulombPeriodic*)this)->ewald = std::make_shared<EwaldPeriodic>(gInfo, pointCharges.size());
	return ewald->energyAndGrad(pointCharges);
}
