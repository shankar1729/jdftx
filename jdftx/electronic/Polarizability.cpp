/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kathleen Schwarz

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

#include <electronic/Polarizability.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>

Polarizability::Polarizability() : eigenBasis(NonInteracting), Ecut(0), nEigs(0)
{
}

inline void pairDensity_thread(int bStart, int bStop, int nV, int nC, const ColumnBundle* C, const diagMatrix* eig, ColumnBundle* rho, diagMatrix* eigDiff)
{	int b = bStart;
	int v = b / nC;
	int c = b - v*nC;
	complexDataRptr conjICv = conj(I(C->getColumn(v)));
	while(b<bStop)
	{	rho->setColumn(b, J(conjICv * I(C->getColumn(nV+c))));
		eigDiff->at(b) = eig->at(nV+c) - eig->at(v);
		//Next cv pair:
		b++; if(b==bStop) break;
		c++;
		if(c==nC) { c=0; v++; conjICv = conj(I(C->getColumn(v))); }
	}
}

inline void coulomb_thread(int bStart, int bStop, const Everything* e, vector3<> dk, const ColumnBundle* rho, ColumnBundle* Krho)
{	for(int b=bStart; b<bStop; b++)
		Krho->setColumn(b, (*(e->coulomb))(rho->getColumn(b), dk, 0.));
}

inline void exCorr_thread(int bStart, int bStop, const DataRptr* exc_nn, const ColumnBundle* rho, ColumnBundle* KXCrho)
{	for(int b=bStart; b<bStop; b++)
		KXCrho->setColumn(b, J((*exc_nn) * I(rho->getColumn(b))));
}

inline void ldaDblPrime_thread(int iStart, int iStop, const double* n, double* e_nn)
{	for(int i=iStart; i<iStop; i++)
	{	double rs = pow((4.*M_PI/3) * n[i], -1./3);
		double rsSqrt = sqrt(rs);
		double rs2 = pow(rs, 2);
		if(rs > 1.)
			e_nn[i] = -0.06105711778985891 * rs2
				* (1.2067926731475986 + rsSqrt)
				* (1.7564571615819247 + rsSqrt)
				* (3.6929643467743336 + 1.5452156337532634*rsSqrt + rs)
				* (1.78462362335091 + 2.395514389326805*rsSqrt + rs)
				/ pow(1. + 1.0529*rsSqrt + 0.3334*rs, 3);
		else
			e_nn[i] = rs2 *
				( 0.009866928037941276 * (-11.754965794910463 + rs) * (7.354022398684049 + rs)
				- 0.001861684535460618 * rs2 * log(rs) );
	}
}

void Polarizability::dump(const Everything& e)
{	
	logPrintf("Dumping polarizability matrix:\n"); logFlush();
	
	if(e.eInfo.qnums.size()>1) die("\nPolarizability currently implemented only for single k-point calculations");
	if(e.exCorr.getName() != "lda-PZ") die("\nPolarizability currently implemented only for lda-PZ exchange correlation");
	if(dk.length()) die("\nPolarizability currently implemented only for dk = 0\n");
	
	if(Ecut<=0.) Ecut = 4.*e.cntrl.Ecut;
	logPrintf("\tSetting up reduced basis at Ecut=%lg: ", Ecut);
	Basis basis; basis.setup(e.gInfo, e.iInfo, Ecut, dk);
	
	logPrintf("\tComputing occupied x unoccupied (CV) pair-densities\n"); logFlush();
	int nV = e.eInfo.nElectrons/2;
	int nC = e.eInfo.nBands - nV;
	if(nC <= 0) die("\nNo unoccupied states available for polarizability calculation.\n");
	int nCV = nC * nV;
	ColumnBundle rho(nCV, basis.nbasis, &basis); //pair density
	diagMatrix eigDiff(nCV); //eigen-value differences (C-V)
	threadOperators = false;
	threadLaunch(pairDensity_thread, nCV, nV, nC, &e.eVars.C[0], &e.eVars.Hsub_eigs[0], &rho, &eigDiff);
	threadOperators = true;
	
	logPrintf("\tApplying Coulomb kernel\n"); logFlush();
	matrix K;
	{	ColumnBundle Krho = rho.similar();
		threadOperators = false;
		threadLaunch(coulomb_thread, nCV, &e, dk, &rho, &Krho);
		threadOperators = true;
		logPrintf("\tForming Coulomb matrix in CV basis\n"); logFlush();
		K = e.gInfo.detR * (rho^Krho);
	}

	logPrintf("\tApplying Exchange-Correlation kernel\n"); logFlush();
	matrix KXC;
	{	DataRptr exc_nn; nullToZero(exc_nn, e.gInfo);
		threadLaunch(ldaDblPrime_thread, e.gInfo.nr, e.eVars.get_nTot()->data(), exc_nn->data());
	
		ColumnBundle KXCrho = rho.similar();
		threadOperators = false;
		threadLaunch(exCorr_thread, nCV, &exc_nn, &rho, &KXCrho);
		threadOperators = true;
		logPrintf("\tForming Exchange-Correlation matrix in CV basis\n"); logFlush();
		KXC = e.gInfo.detR * (rho^KXCrho);
	}
	
	//Compute operator matrices in current (CV) basis
	logPrintf("\tComputing polarizability matrices in CV basis\n"); logFlush();
	matrix invXni = -0.25*eigDiff; //inverse of non-interacting susceptibility
	matrix invXtot = invXni - KXC; //inverse of charge response to total electrostatic potential
	matrix invXext = invXtot - K; //inverse of charge response to external electrostatic potential

	//Determine transformation to chosen eigen-basis
	logPrintf("\tComputing transformation matrix from CV to chosen eigen-basis\n"); logFlush();
	if(nEigs<=0) nEigs = nCV;
	matrix Q; //transformation matrix from CV to chosen basis
	{	matrix Umhalf = invsqrt(e.gInfo.detR*(rho^rho)); //matrix that orthogonalizes CV basis
		const matrix* invXbasis = 0;
		switch(eigenBasis)
		{	case NonInteracting: invXbasis = &invXni; break;
			case External: invXbasis = &invXext; break;
			case Total: invXbasis = &invXtot; break;
			default: assert(!"Invalid eigenBasis");
		}
		matrix evecs; diagMatrix eigs;
		(Umhalf * (*invXbasis) * Umhalf).diagonalize(evecs, eigs);
		//Truncate eigen-expansion:
		//--- most negative eigenvalue of neg-definite Xbasis is least negative eigenvalue of invXbasis
		//--- eigenvalues in ascending order => pick last nEigs eigenvalues
		Q = Umhalf * evecs(0,nCV, nCV-nEigs,nCV);
	}

	//Transform all quantities to eigenbasis:
	logPrintf("\tTransforming output quantities to chosen eigen-basis\n"); logFlush();
	rho = rho * Q; //now equal to what we call V in derivations
	invXni = dagger(Q) * invXni * Q;
	invXext = dagger(Q) * invXext * Q;
	invXtot = dagger(Q) * invXtot * Q;
	K = dagger(Q) * K * Q;
	KXC = dagger(Q) * KXC * Q;

	//Dump:
	logPrintf("\tDumping '%s' ... ", e.dump.getFilename("pol_*").c_str()); logFlush();
	rho.write(e.dump.getFilename("pol_basis").c_str());
	invXni.write(e.dump.getFilename("pol_invXni").c_str());
	invXext.write(e.dump.getFilename("pol_invXext").c_str());
	invXtot.write(e.dump.getFilename("pol_invXtot").c_str());
	K.write(e.dump.getFilename("pol_K").c_str());
	KXC.write(e.dump.getFilename("pol_KXC").c_str());
	logPrintf("Done.\n");
	logFlush();
}
