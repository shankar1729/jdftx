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
#include <core/LatticeUtils.h>
#include <core/DataIO.h>

Polarizability::Polarizability() : eigenBasis(NonInteracting), Ecut(0), nEigs(0)
{
}


class PairDensityCalculator
{
	const Everything& e; int ik, nK;
	const ColumnBundle *C1, *C2;
	const diagMatrix *eig1, *eig2;
	int *index1, *index2;
	
	//Setup transformed index array in targetIndex
	void setupIndex(const Supercell::KmeshTransform kTransform, const Basis& basis, int* &targetIndex) const
	{	const matrix3<int> mRot = (~e.symm.getMatrices()[kTransform.iSym]) * kTransform.invert;
		std::vector<int> index(basis.nbasis);
		for(unsigned j=0; j<basis.nbasis; j++)
			index[j] = e.gInfo.fullGindex(mRot * basis.iGarr[j] - kTransform.offset);
		#ifdef GPU_ENABLED
		cudaMalloc(&targetIndex, sizeof(int)*index.size());
		cudaMemcpy(targetIndex, index.data(), sizeof(int)*index.size(), cudaMemcpyHostToDevice);
		#else
		targetIndex = new int[index.size()];
		eblas_copy(targetIndex, index.data(), index.size());
		#endif
	}
	
	//ColumnBundle::getColumn with custom index array:
	static complexDataGptr getColumn(const ColumnBundle* C, int col, const int* index)
	{	complexDataGptr full; nullToZero(full, *(C->basis->gInfo)); //initialize a full G-space vector to zero
		callPref(eblas_scatter_zdaxpy)(C->basis->nbasis, 1., index, C->dataPref()+C->index(col,0), full->dataPref()); //scatter from col'th column
		return full;
	}
	
public:
	//Setup to compute pair densities between kmesh[ik] and its partner dk away
	PairDensityCalculator(const Everything& e, const vector3<>& dk, int ik) : e(e), ik(ik)
	{
		//Find the transformations for the two k-points:
		const std::vector< vector3<> >& kmesh = e.coulombParams.supercell->kmesh;
		const std::vector<Supercell::KmeshTransform>& kmeshTransform = e.coulombParams.supercell->kmeshTransform;
		Supercell::KmeshTransform kTransform1 = kmeshTransform[ik], kTransform2;
		nK = kmesh.size();
		vector3<> k2 = kmesh[ik] + dk;
		bool foundk2 = false;
		for(unsigned ik2=0; ik2<kmesh.size(); ik2++)
			if(circDistanceSquared(kmesh[ik2],k2) < symmThresholdSq)
			{	kTransform2 = kmeshTransform[ik2];
				vector3<> extraOffsetTemp = k2 - kmesh[ik2];
				vector3<int> extraOffset;
				for(int k=0; k<3; k++)
				{	extraOffset[k] = int(round(extraOffsetTemp[k]));
					assert(fabs(extraOffset[k]-extraOffsetTemp[k]) < symmThreshold);
				}
				kTransform2.offset += extraOffset;
				foundk2 = true;
				break;
			}
		assert(foundk2); //such a partner should always be found for a uniform kmesh
		
		//Get pointers to source wavefunctons and eigenvalues:
		C1 = &(e.eVars.C[kTransform1.iReduced]);
		C2 = &(e.eVars.C[kTransform2.iReduced]);
		eig1 = &(e.eVars.Hsub_eigs[kTransform1.iReduced]);
		eig2 = &(e.eVars.Hsub_eigs[kTransform2.iReduced]);
		
		//Setup index arrays for the two k-points:
		setupIndex(kTransform1, *(C1->basis), index1);
		setupIndex(kTransform2, *(C2->basis), index2);
	}
	
	~PairDensityCalculator()
	{
		#ifdef GPU_ENABLED
		cudaFree(index1); cudaFree(index2);
		#else
		delete[] index1; delete[] index2;
		#endif
	}

	//Store resulting pair densities scaled by 2*invsqrt(eigenvalue differences) in rho,
	//so that the non-interacting susceptibility is negative identity in this basis.
	void compute(int nV, int nC, ColumnBundle& rho, int kOffset) const
	{	threadLaunch(isGpuEnabled() ? 1 : 0, compute_thread, nV*nC, nV, nC, &rho, kOffset, this);
	}
	
	//Accumulate contribution from currentkpoint pair to negative of noninteracting susceptibility in plane-wave basis:
	void accumMinusXniPW(int nV, int nC, const Basis& basis, matrix& minusXni)
	{	assert(minusXni.nRows() == int(basis.nbasis));
		assert(minusXni.nCols() == int(basis.nbasis));
		ColumnBundle rho(nV*nC, basis.nbasis, &basis);
		compute(nV, nC, rho, 0);
		//Xni += (detR)*rho*dagger(rho):
		callPref(eblas_zgemm)(CblasNoTrans, CblasConjTrans, basis.nbasis, basis.nbasis, rho.nCols(),
			basis.gInfo->detR, rho.dataPref(), rho.colLength(), rho.dataPref(), rho.colLength(),
			1., minusXni.dataPref(), minusXni.nRows());
	}
	
private:
	void compute_sub(int bStart, int bStop, int nV, int nC, ColumnBundle* rho, int kOffset) const
	{	int b = bStart;
		int v = b / nC;
		int c = b - v*nC;
		complexDataRptr conjICv = conj(I(getColumn(C1, v, index1)));
		while(b<bStop)
		{	double sqrtEigDen = sqrt(4./(nK * (eig2->at(nV+c) - eig1->at(v))));
			rho->setColumn(kOffset+b, sqrtEigDen * J(conjICv * I(getColumn(C2, nV+c, index2))));
			//Next cv pair:
			b++; if(b==bStop) break;
			c++;
			if(c==nC) { c=0; v++; conjICv = conj(I(getColumn(C1, v, index1))); }
		}
	}
	static void compute_thread(int bStart, int bStop, int nV, int nC, ColumnBundle* rho, int kOffset, const PairDensityCalculator* pdc)
	{	pdc->compute_sub(bStart, bStop, nV, nC, rho, kOffset);
	}
};

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
	{	if(n[i]<1e-4) { e_nn[i] = 0.; continue; }
		double rs = pow((4.*M_PI/3) * n[i], -1./3);
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
	if(e.exCorr.getName() != "lda-PZ") die("\nPolarizability currently implemented only for lda-PZ exchange correlation");
	
	if(Ecut<=0.) Ecut = 4.*e.cntrl.Ecut;
	logPrintf("\tSetting up reduced basis at Ecut=%lg: ", Ecut);
	Basis basis; basis.setup(e.gInfo, e.iInfo, Ecut, dk);
	
	int nV = e.eInfo.nElectrons/2;
	int nC = e.eInfo.nBands - nV;
	int nK = e.coulombParams.supercell->kmesh.size();
	if(nC <= 0) die("\nNo unoccupied states available for polarizability calculation.\n");
	int nCVK = nC * nV * nK;
	
	//Determine whether to start out in CV or PW basis:
	bool pwBasis = (2*nCVK > int(basis.nbasis)); //switch to PW basis a little early since CV basis begins to become numerically unstable
	int nColumns = pwBasis ? int(basis.nbasis) : nCVK;
	const char* basisName = pwBasis ? "PW" : "CV";
	
	ColumnBundle V(nColumns, basis.nbasis, &basis); //orthonormal basis vectors
	matrix Xni; //non-interacting susceptibility (in basis V)
	
	if(pwBasis)
	{	logPrintf("\tComputing NonInteracting polarizability in plane-wave basis\n"); logFlush();
		//Set the basis to an identity matrix:
		V.zero();
		complex* Vdata = V.data();
		double invsqrtVol = 1./sqrt(e.gInfo.detR);
		for(int i=0; i<nColumns; i++) Vdata[V.index(i,i)] = invsqrtVol;
		//Get the PW basis non-interacting susceptibility matrix:
		matrix minusXni(nColumns, nColumns); minusXni.zero();
		threadOperators = false;
		for(int ik=0; ik<nK; ik++)
			PairDensityCalculator(e, dk, ik).accumMinusXniPW(nV, nC, basis, minusXni);
		threadOperators = true;
		Xni = -minusXni;
	}
	else
	{	logPrintf("\tComputing occupied x unoccupied (CV) pair-densities and NonInteracting polarizability\n"); logFlush();
		diagMatrix eigDiff(nColumns); //eigen-value differences (C-V)
		threadOperators = false;
		for(int ik=0; ik<nK; ik++)
			PairDensityCalculator(e, dk, ik).compute(nV, nC, V, ik*nV*nC);
		threadOperators = true;
		matrix invXni = -eye(nColumns); //inverse of non-interacting susceptibility
		logPrintf("\tOrthonormalizing basis\n"); logFlush();
		matrix Umhalf = invsqrt(e.gInfo.detR*(V^V));
		V = V * Umhalf;
		Xni = dagger_symmetrize(inv(Umhalf * invXni * Umhalf));
	}
	
	logPrintf("\tApplying Coulomb kernel\n"); logFlush();
	matrix K;
	{	ColumnBundle KV = V.similar();
		threadOperators = false;
		threadLaunch(isGpuEnabled() ? 1 : 0, coulomb_thread, nColumns, &e, dk, &V, &KV);
		threadOperators = true;
		logPrintf("\tForming Coulomb matrix in %s basis\n", basisName); logFlush();
		K = e.gInfo.detR * (V^KV);
		K.write(e.dump.getFilename("pol_K").c_str());
	}

	logPrintf("\tApplying Exchange-Correlation kernel\n"); logFlush();
	matrix KXC;
	{	DataRptr exc_nn; nullToZero(exc_nn, e.gInfo);
		threadLaunch(ldaDblPrime_thread, e.gInfo.nr, e.eVars.get_nTot()->data(), exc_nn->data());
		saveRawBinary(exc_nn, e.dump.getFilename("pol_ExcDblPrime").c_str());
		
		ColumnBundle KXCV = V.similar();
		threadOperators = false;
		threadLaunch(isGpuEnabled() ? 1 : 0, exCorr_thread, nColumns, &exc_nn, &V, &KXCV);
		threadOperators = true;
		logPrintf("\tForming Exchange-Correlation matrix in %s basis\n", basisName); logFlush();
		KXC = e.gInfo.detR * (V^KXCV);
	}
	
	//Compute operator matrices in current (CV) basis
	logPrintf("\tComputing External and Total polarizability matrices in %s basis\n", basisName); logFlush();
	matrix Xtot = dagger_symmetrize(inv(eye(nColumns) - Xni*(  KXC  )) * Xni); //charge response to total electrostatic potential
	matrix Xext = dagger_symmetrize(inv(eye(nColumns) - Xni*(K + KXC)) * Xni); //charge response to external electrostatic potential
	
	//Compute dielectric band structure:
	{	string fname = e.dump.getFilename("epsInvEigs");
		logPrintf("\tDumping '%s' ... ", fname.c_str()); logFlush();
		matrix epsInvEvecs; diagMatrix epsInvEigs;
		matrix Khalf = pow(dagger_symmetrize(K), 0.5);
		(eye(nColumns) + Khalf * Xext * Khalf).diagonalize(epsInvEvecs, epsInvEigs); //epsInv (symmetrized)
		FILE* fp=fopen(fname.c_str(), "w");
		epsInvEigs.print(fp, "%.15f\n");
		fclose(fp);
		logPrintf("Done.\n"); logFlush();
	}

	//Determine transformation to chosen eigen-basis
	extern EnumStringMap<EigenBasis> polarizabilityMap;
	logPrintf("\tComputing transformation matrix from %s to %s polarizability eigen-basis\n", basisName, polarizabilityMap.getString(eigenBasis)); logFlush();
	const matrix* Xbasis = 0;
	switch(eigenBasis)
	{	case NonInteracting: Xbasis = &Xni; break;
		case External: Xbasis = &Xext; break;
		case Total: Xbasis = &Xtot; break;
		default: assert(!"Invalid eigenBasis");
	}
	if(nEigs<=0 || nEigs>nColumns) nEigs = nColumns;
	matrix Q; //transformation matrix from CV to chosen basis
	{	matrix evecs; diagMatrix eigs;
		(*Xbasis).diagonalize(evecs, eigs);
		Q = evecs(0,nColumns, 0,nEigs); //Largest negative eigenvalues of Xbasis appear in the beginning; select first nEigs of them
	}
	
	//Transform all quantities to eigenbasis:
	logPrintf("\tTransforming output quantities to %s polarizability eigen-basis\n", polarizabilityMap.getString(eigenBasis)); logFlush();
	V = V * Q;
	Xni = dagger(Q) * Xni * Q;
	Xext = dagger(Q) * Xext * Q;
	Xtot = dagger(Q) * Xtot * Q;
	K = dagger(Q) * K * Q;
	KXC = dagger(Q) * KXC * Q;

	//Dump:
	logPrintf("\tDumping '%s' ... ", e.dump.getFilename("pol_*").c_str()); logFlush();
	V.write(e.dump.getFilename("pol_basis").c_str());
	Xni.write(e.dump.getFilename("pol_Xni").c_str());
	Xext.write(e.dump.getFilename("pol_Xext").c_str());
	Xtot.write(e.dump.getFilename("pol_Xtot").c_str());
	K.write(e.dump.getFilename("pol_K").c_str());
	KXC.write(e.dump.getFilename("pol_KXC").c_str());
	logPrintf("Done.\n");
	logFlush();
}
