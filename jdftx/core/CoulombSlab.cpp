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

#include <core/CoulombSlab.h>
#include <core/Coulomb_internal.h>
#include <core/BlasExtra.h>

//! 2D Ewald sum
struct EwaldSlab
{
	const GridInfo& gInfo;
	int iDir; //!< truncated direction
	double ionMargin; //!< margin between min image distance between ions and half truncation length
	double sigma; //!< gaussian width for Ewald sums
	vector3<int> Nreal; //!< max unit cell indices for real-space sum
	vector3<int> Nrecip; //!< max unit cell indices for reciprocal-space sum
	
	EwaldSlab(const GridInfo& gInfo, int iDir, double ionMargin) : gInfo(gInfo), iDir(iDir), ionMargin(ionMargin)
	{	logPrintf("\n---------- Setting up 2D ewald sum ----------\n");
		//Determine optimum gaussian width for 2D Ewald sums:
		// From below, the number of reciprocal cells ~ Prod_{k!=iDir} |R.column[k]|
		//    and number of real space cells ~ Prod_{k!=iDir} |G.row[k]|
		// including the fact that the reciprocal space sum has 5 times
		// as many special function calls per term as the real space sum
		sigma = 1.;
		for(int k=0; k<3; k++)
			if(k != iDir)
				sigma *= gInfo.R.column(k).length() / gInfo.G.row(k).length();
		sigma = pow(5.*sigma, 1./4);
		logPrintf("Optimum gaussian width for ewald sums = %lf bohr.\n", sigma);
		
		//Carry real space sums to Rmax = 10 sigma and Gmax = 10/sigma
		//This leads to relative errors ~ 1e-22 in both sums, well within double precision limits
		for(int k=0; k<3; k++)
		{	Nreal[k]  = (k==iDir) ? 0 : 1+ceil(10. * gInfo.G.row(k).length() * sigma / (2*M_PI));
			Nrecip[k] = (k==iDir) ? 0 : 1+ceil(10. * gInfo.R.column(k).length() / (2*M_PI*sigma));
		}
		logPrintf("Real space sums over %d unit cells with max indices ", (2*Nreal[0]+1)*(2*Nreal[1]+1)*(2*Nreal[2]+1));
		Nreal.print(globalLog, " %d ");
		logPrintf("Reciprocal space sums over %d unit cells with max indices ", (2*Nrecip[0]+1)*(2*Nrecip[1]+1)*(2*Nrecip[2]+1));
		Nrecip.print(globalLog, " %d ");
	}
	
	double energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
	{	if(!pointCharges.size()) return 0.;
		double eta = sqrt(0.5)/sigma, etaSq=eta*eta, etaSqrtPiInv = 1./(eta*sqrt(M_PI));
		double sigmaSq = sigma * sigma;
		//Position independent terms: (Self-energy correction)
		double ZsqTot = 0.;
		for(const Coulomb::PointCharge& pc: pointCharges)
			ZsqTot += pc.Z * pc.Z;
		double E = -0.5 * ZsqTot * eta * (2./sqrt(M_PI));
		//Reduce positions to first unit cell:
		//Shift all points in the truncated direction into the 1D Wigner-Seitz cell
		//centered on one of the atoms; choice of this atom is irrelevant if every atom
		//lies in the WS cell of the other with a consistent translation:
		vector3<> pos0(0.,0.,0.);
		pos0[iDir] = pointCharges[0].pos[iDir];
		for(Coulomb::PointCharge& pc: pointCharges)
			for(int k=0; k<3; k++)
				pc.pos[k] -= floor(0.5 + pc.pos[k] - pos0[k]);
		//Real space sum:
		vector3<int> iR; //integer cell number
		for(const Coulomb::PointCharge& pc2: pointCharges)
			for(Coulomb::PointCharge& pc1: pointCharges)
				for(iR[0]=-Nreal[0]; iR[0]<=Nreal[0]; iR[0]++)
					for(iR[1]=-Nreal[1]; iR[1]<=Nreal[1]; iR[1]++)
						for(iR[2]=-Nreal[2]; iR[2]<=Nreal[2]; iR[2]++)
						{	vector3<> x = iR + (pc1.pos - pc2.pos);
							double rSq = gInfo.RTR.metric_length_squared(x);
							if(!rSq) continue; //exclude self-interaction
							double r = sqrt(rSq);
							E += 0.5 * pc1.Z * pc2.Z * erfc(eta*r)/r;
							pc1.force += (gInfo.RTR * x) *
								(pc1.Z * pc2.Z * (erfc(eta*r)/r + (2./sqrt(M_PI))*eta*exp(-etaSq*rSq))/rSq);
						}
		//Reciprocal space sum:
		double L = sqrt(gInfo.RTR(iDir,iDir)); //length of truncated direction
		double volPrefac = M_PI * L / gInfo.detR;
		for(unsigned i1=0; i1<pointCharges.size(); i1++)
		{	Coulomb::PointCharge& pc1 = pointCharges[i1];
			for(unsigned i2=0; i2<=i1; i2++)
			{	Coulomb::PointCharge& pc2 = pointCharges[i2];
				double prefac = volPrefac * pc1.Z * pc2.Z * (i1==i2 ? 1 : 2);
				vector3<> r12 = pc1.pos - pc2.pos;
				double z12 = L * r12[iDir];
				if(fabs(z12) >= 0.5*L-ionMargin)
					die("Separation between atoms %d and %d lies in the truncation margin.\n", i1, i2);
				double E12 = 0.; vector3<> E12_r12(0.,0.,0.); //energy and gradient from this pair
				vector3<int> iG; //integer reciprocal cell number (iG[iDir] will remain 0)
				for(iG[0]=-Nrecip[0]; iG[0]<=Nrecip[0]; iG[0]++)
					for(iG[1]=-Nrecip[1]; iG[1]<=Nrecip[1]; iG[1]++)
						for(iG[2]=-Nrecip[2]; iG[2]<=Nrecip[2]; iG[2]++)
						{	//2D structure factor term and derivative
							double c, s; sincos((2*M_PI)*dot(iG,r12), &s, &c);
							//Contribution from truncated direction:
							double Gsq = gInfo.GGT.metric_length_squared(iG);
							double zTerm, zTermPrime;
							if(Gsq)
							{	double G = sqrt(Gsq);
								if(fabs(G*z12) > 100.) continue; //negligible contribution
								double expPlus = exp(G*z12), expMinus = 1./expPlus;
								double erfcPlus  = erfc(eta*(sigmaSq*G + z12));
								double erfcMinus = erfc(eta*(sigmaSq*G - z12));
								zTerm = (0.5/G) * (expPlus * erfcPlus + expMinus * erfcMinus);
								zTermPrime = 0.5 * (expPlus * erfcPlus - expMinus * erfcMinus);
							}
							else
							{	double erfz = erf(eta * z12);
								double gauss = exp(-etaSq * z12*z12);
								zTerm = -z12 * erfz - etaSqrtPiInv * gauss;
								zTermPrime = -erfz;
							}
							//Update energy and forces:
							E12 += prefac * c * zTerm;
							E12_r12 += (prefac * -s * zTerm * (2*M_PI)) * iG;
							E12_r12[iDir] += prefac * c * zTermPrime * L;
						}
				E += E12;
				pc1.force -= E12_r12;
				pc2.force += E12_r12;
			}
		}
		return E;
	}
};


CoulombSlab::CoulombSlab(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params)
{	int iDir = params.iDir;
	string dirName(3, '0'); dirName[iDir] = '1';
	if(fabs(gInfo.GGT(iDir,iDir) * gInfo.RTR(iDir,iDir) - 4*M_PI*M_PI) > 1e-14)
		die("Lattice direction %s is not perpendicular to the other two basis vectors.\n", dirName.c_str());
	logPrintf("Initialized slab truncation along lattice direction %s\n", dirName.c_str());
}

DataGptr CoulombSlab::operator()(DataGptr&& in) const
{	int iDir = params.iDir;
	double hlfL = 0.5*sqrt(gInfo.RTR(iDir,iDir));
	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombSlab_calc(iDir, hlfL), in->dataPref(false));
	return in;
}

double CoulombSlab::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	if(!ewald)
		((CoulombSlab*)this)->ewald = std::make_shared<EwaldSlab>(gInfo, params.iDir, params.ionMargin);
	return ewald->energyAndGrad(pointCharges);
}
