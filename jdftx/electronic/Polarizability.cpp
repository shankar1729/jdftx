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
#include <electronic/ColumnBundleTransform.h>
#include <core/LatticeUtils.h>
#include <core/VectorField.h>
#include <core/ScalarFieldIO.h>

Polarizability::Polarizability() : eigenBasis(NonInteracting), Ecut(0), nEigs(0)
{
}


class PairDensityCalculator
{
	int nK;
	
	struct State
	{	const ColumnBundle* C;
		const diagMatrix* eig;
		std::shared_ptr<ColumnBundleTransform> transform;
		
		void setup(const Everything& e, vector3<> k, const Supercell::KmeshTransform& kTransform)
		{	//Get the columnbundle and eigenvalues:
			C = &(e.eVars.C[kTransform.iReduced]);
			eig = &(e.eVars.Hsub_eigs[kTransform.iReduced]);
			//Compute the index array
			logSuspend();
			basisOut.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, k);
			logResume();
			transform = std::make_shared<ColumnBundleTransform>(C->qnum->k, *(C->basis), k, basisOut,
				e.eInfo.spinorLength(), e.symm.getMatrices()[kTransform.iSym], kTransform.invert);
		}
		
		void setup(const Everything& e, vector3<> k, string fnameWfns, string fnameEig)
		{	//Wrap k point to (-0.5,0.5] (consistent with ElecInfo):
			vector3<> kWrapped = k;
			for(int j=0; j<3; j++)
				kWrapped[j] -= ceil(kWrapped[j]-0.5);
			//Setup basis:
			logSuspend();
			basisExt.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, kWrapped);
			logResume();
			//Read wavefunctions:
			Cext.init(e.eInfo.nBands, basisExt.nbasis*e.eInfo.spinorLength(), &basisExt, 0);
			off_t flen = fileSize(fnameWfns.c_str());
			int nBytesPerBand = basisExt.nbasis*sizeof(complex);
			if(flen < long(Cext.nData()*sizeof(complex))) die("\nFile '%s' does not exist or is too short.\n", fnameWfns.c_str());
			if(flen % nBytesPerBand) die("\nFile '%s' is not a multiple of %d bytes per band (basis mismatch?).\n", fnameWfns.c_str(), nBytesPerBand);
			((ManagedMemory<complex>&)Cext).read(fnameWfns.c_str());
			C = &Cext;
			//Read eigenvalues:
			eigExt.resize(e.eInfo.nBands);
			if(fileSize(fnameEig.c_str()) < 0) die("\nFile '%s' does not exist.\n", fnameEig.c_str());
			FILE* fpEig = fopen(fnameEig.c_str(), "r");
			freadLE(eigExt.data(), sizeof(double), eigExt.nRows(), fpEig);
			if(feof(fpEig))  die("\nFile '%s' ended before all eigenvalues could be read.\n", fnameEig.c_str());
			fclose(fpEig);
			eig = &eigExt;
			//Setup the output basis and transformation:
			logSuspend();
			basisOut.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, k);
			logResume();
			transform = std::make_shared<ColumnBundleTransform>(kWrapped, basisExt, k, basisOut,
				e.eInfo.spinorLength(), SpaceGroupOp(), +1);
		}
		
		//ColumnBundle::getColumn, but with custom index array:
		complexScalarFieldTilde getColumn(int col) const
		{	ColumnBundle Cout(1, C->colLength(), &basisOut, 0, isGpuEnabled());
			Cout.zero();
			transform->scatterAxpy(1., *C,col, Cout,0);
			return Cout.getColumn(0,0);
		}
		
	private:
		ColumnBundle Cext; diagMatrix eigExt; Basis basisExt, basisOut; //Externally read-in wavefunction, corresponding basis and eigenvalues
	}
	state1, state2;

public:
	//Setup to compute pair densities between kmesh[ik] and its partner dk away
	PairDensityCalculator(const Everything& e, const vector3<>& dk, int ik)
	{
		//Find the transformations / data sources for the two k-points:
		const std::vector< vector3<> >& kmesh = e.coulombParams.supercell->kmesh;
		const std::vector<Supercell::KmeshTransform>& kmeshTransform = e.coulombParams.supercell->kmeshTransform;
		nK = kmesh.size();
		vector3<> k2 = kmesh[ik] + dk;
		
		state1.setup(e, kmesh[ik], kmeshTransform[ik]); //setup first state (always from current system's kmesh)
		
		if(e.dump.polarizability->dkFilenamePattern.length()) //get second state from external data
		{	state2.setup(e, k2,
				e.dump.polarizability->dkFilename(ik,"wfns"),
				e.dump.polarizability->dkFilename(ik,"eigenvals") );
		}
		else //get second state from current system's kmesh as well
		{	bool foundk2 = false;
			Supercell::KmeshTransform kTransform2;
			for(unsigned ik2=0; ik2<kmesh.size(); ik2++)
				if(circDistanceSquared(kmesh[ik2],k2) < symmThresholdSq)
				{	double offsetErr;
					kTransform2 = kmeshTransform[ik2];
					kTransform2.offset += round(k2 - kmesh[ik2], &offsetErr);
					foundk2 = true;
					assert(offsetErr < symmThreshold);
					break;
				}
			assert(foundk2); //such a partner should always be found for a uniform kmesh
			state2.setup(e, k2, kTransform2);
		}
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
		complexScalarField conjICv = conj(I(state1.getColumn(v)));
		while(b<bStop)
		{	double sqrtEigDen = sqrt(4./(nK * (state2.eig->at(nV+c) - state1.eig->at(v))));
			rho->setColumn(kOffset+b,0, sqrtEigDen * J(conjICv * I(state2.getColumn(nV+c))));
			//Next cv pair:
			b++; if(b==bStop) break;
			c++;
			if(c==nC) { c=0; v++; conjICv = conj(I(state1.getColumn(v))); }
		}
	}
	static void compute_thread(int bStart, int bStop, int nV, int nC, ColumnBundle* rho, int kOffset, const PairDensityCalculator* pdc)
	{	pdc->compute_sub(bStart, bStop, nV, nC, rho, kOffset);
	}
};

inline void coulomb_thread(int bStart, int bStop, const Everything* e, vector3<> dk, const ColumnBundle* rho, ColumnBundle* Krho)
{	const GridInfo& gInfoWfns = *(rho->basis->gInfo);
	for(int b=bStart; b<bStop; b++)
	{	complexScalarFieldTilde rho_b = rho->getColumn(b,0);
		if(&gInfoWfns != &e->gInfo) rho_b = changeGrid(rho_b, e->gInfo);
		complexScalarFieldTilde Krho_b = (*(e->coulomb))(rho_b, dk, 0.);
		if(&gInfoWfns != &e->gInfo) Krho_b = changeGrid(Krho_b, gInfoWfns);
		Krho->setColumn(b,0, Krho_b);
	}
}

matrix coulombMatrix(const ColumnBundle& V, const Everything& e, vector3<> dk)
{	ColumnBundle KV = V.similar();
	threadLaunch(isGpuEnabled() ? 1 : 0, coulomb_thread, V.nCols(), &e, dk, &V, &KV);
	logPrintf("\tForming Coulomb matrix\n"); logFlush();
	return e.gInfo.detR * (V^KV);
}


//------- Exchange and correlation -----------
typedef ScalarFieldMultiplet<complexScalarFieldData,3> complexScalarFieldVec;

//Get the gradient of one column of a column bundle
complexScalarFieldVec gradient(const ColumnBundle& Y, int col)
{	ColumnBundle Ysub = Y.getSub(col, col+1);
	complexScalarFieldVec DY;
	for(int j=0; j<3; j++)
		DY[j] = I(D(Ysub,j).getColumn(0,0));
	return DY;
}

//Accumulate the divergence of a complex vector field into one column of a columnbundle
void axpyDivergence(double alpha, const complexScalarFieldVec& x, ColumnBundle& Y, int col)
{	ColumnBundle Ysub = Y.getSub(col, col+1);
	ColumnBundle xj = Ysub.similar();
	for(int j=0; j<3; j++)
	{	xj.setColumn(0,0, J(x[j]));
		Ysub += alpha * D(xj, j);
	}
	Y.setSub(col, Ysub);
}

complexScalarField dotElemwise(const VectorField& x, const complexScalarFieldVec& y)
{	complexScalarField ret;
	for(int j=0; j<3; j++) ret += x[j] * y[j];
	return ret;
}

complexScalarFieldVec operator*(const VectorField& x, const complexScalarField& y)
{	complexScalarFieldVec ret;
	for(int j=0; j<3; j++) ret[j] = x[j] * y;
	return ret;
}

complexScalarFieldVec operator*(const ScalarField& x, const complexScalarFieldVec& y)
{	complexScalarFieldVec ret;
	for(int j=0; j<3; j++) ret[j] = x * y[j];
	return ret;
}

inline void exCorr_thread(int bStart, int bStop, const ScalarField* exc_nn, const VectorField* Dn,
	const ScalarField* exc_sigma, const ScalarField* exc_nsigma, const ScalarField* exc_sigmasigma,
	const ColumnBundle* rho, ColumnBundle* KXCrho)
{	for(int b=bStart; b<bStop; b++)
	{	//Get the basis vector (and optionally its gradient) in real space:
		complexScalarField V = I(rho->getColumn(b,0)); complexScalarFieldVec DV;
		if(*exc_sigma) DV = gradient(*rho, b);
		//Add contributions which are local towards the right:
		complexScalarField KV = (*exc_nn) * V, DnDV;
		if(*exc_sigma)
		{	DnDV = 2. * dotElemwise(*Dn, DV);
			KV += (*exc_nsigma) * DnDV;
		}
		KXCrho->setColumn(b,0, J(KV));
		//Add contributions which have a gradient towards the right:
		if(*exc_sigma)
		{	complexScalarField DnTerm = (*exc_nsigma)*V + (*exc_sigmasigma)*DnDV;
			axpyDivergence(-2., (*Dn) * DnTerm + (*exc_sigma) * DV, *KXCrho, b);
		}
	}
}

matrix exCorrMatrix(const ColumnBundle& V, const Everything& e, const ScalarField& n, vector3<> dk)
{	//Get second derivatives w.r.t density (and gradients)
	ScalarField exc_nn, exc_sigma, exc_nsigma, exc_sigmasigma;
	VectorField Dn;
	e.exCorr.getSecondDerivatives(n, exc_nn, exc_sigma, exc_nsigma, exc_sigmasigma);
	if(exc_sigma) Dn = gradient(n); //needed for GGAs
	//Compute matrix:
	ColumnBundle KXCV = V.similar();
	threadLaunch(isGpuEnabled() ? 1 : 0, exCorr_thread, V.nCols(), &exc_nn, &Dn, &exc_sigma, &exc_nsigma, &exc_sigmasigma, &V, &KXCV);
	logPrintf("\tForming Exchange-Correlation matrix\n"); logFlush();
	return e.gInfo.detR * (V^KXCV);
}


void Polarizability::dump(const Everything& e)
{	
	logPrintf("Dumping polarizability matrix:\n"); logFlush();
	if(e.eInfo.spinType != SpinNone) die("\nPolarizability currently implemented only for spin-unpolarized systems\n");
	
	const std::vector< vector3<> >& kmesh = e.coulombParams.supercell->kmesh;
	if(dkFilenamePattern.length())
	{	//Check if the required files seem to exist - if not print the required info for someone to generate them:
		if(fileSize(dkFilename(0,"wfns").c_str()) <= 0)
		{	logPrintf("\tSave band structure states for the following k-points in the indicated locations:\n");
			for(unsigned ik=0; ik<kmesh.size(); ik++)
			{	vector3<> k2 = kmesh[ik] + dk;
				for(int j=0; j<3; j++) k2[j] -= ceil(k2[j]-0.5); //Reduce to (-0.5,0.5] (consistent with ElecInfo)
				logPrintf("\t\t%20.17f %20.17f %20.17f -> %s\n", k2[0], k2[1], k2[2], dkFilename(ik,"$VAR").c_str());
			}
			logPrintf("\tRerun after generating above files to get polarizability.\n");
			return;
		}
	}
	
	if(Ecut<=0.) Ecut = 4.*e.cntrl.Ecut;
	logPrintf("\tSetting up reduced basis at Ecut=%lg: ", Ecut);
	Basis basis; basis.setup(e.gInfo, e.iInfo, Ecut, dk);
	
	int nV = e.eInfo.nElectrons/2;
	int nC = e.eInfo.nBands - nV;
	int nK = kmesh.size();
	if(nC <= 0) die("\nNo unoccupied states available for polarizability calculation.\n");
	int nCVK = nC * nV * nK;
	
	//Determine whether to start out in CV or PW basis:
	bool pwBasis = (2*nCVK > int(basis.nbasis)); //switch to PW basis a little early since CV basis begins to become numerically unstable
	int nColumns = pwBasis ? int(basis.nbasis) : nCVK;
	const char* basisName = pwBasis ? "PW" : "CV";
	
	QuantumNumber qnum; qnum.k = dk; qnum.spin = 0; qnum.weight = 1./nK;
	ColumnBundle V(nColumns, basis.nbasis, &basis, &qnum); //orthonormal basis vectors
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
		for(int ik=0; ik<nK; ik++)
			PairDensityCalculator(e, dk, ik).accumMinusXniPW(nV, nC, basis, minusXni);
		Xni = -minusXni;
	}
	else
	{	logPrintf("\tComputing occupied x unoccupied (CV) pair-densities and NonInteracting polarizability\n"); logFlush();
		for(int ik=0; ik<nK; ik++)
			PairDensityCalculator(e, dk, ik).compute(nV, nC, V, ik*nV*nC);
		matrix invXni = -eye(nColumns); //inverse of non-interacting susceptibility
		logPrintf("\tOrthonormalizing basis\n"); logFlush();
		matrix Umhalf = invsqrt(e.gInfo.detR*(V^V));
		V = V * Umhalf;
		Xni = dagger_symmetrize(inv(Umhalf * invXni * Umhalf));
	}
	
	logPrintf("\tClearing orthogonal wavefunctions (C) to free memory.\n"); ((Everything&)e).eVars.C.clear();
	
	logPrintf("\tApplying Coulomb kernel\n"); logFlush();
	matrix K = coulombMatrix(V, e, dk);
	
	logPrintf("\tApplying Exchange-Correlation kernel\n"); logFlush();
	matrix KXC = exCorrMatrix(V, e, e.eVars.get_nTot(), dk);
	
	//Compute operator matrices in current (CV) basis
	logPrintf("\tComputing External and Total polarizability matrices in %s basis\n", basisName); logFlush();
	matrix Xtot = dagger_symmetrize(inv(eye(nColumns) - Xni*(  KXC  )) * Xni); //charge response to total electrostatic potential
	matrix Xext = dagger_symmetrize(inv(eye(nColumns) - Xni*(K + KXC)) * Xni); //charge response to external electrostatic potential
	
	//Compute dielectric band structure:
	{	string fname = e.dump.getFilename("epsInvEigs");
		logPrintf("\tDumping '%s' ... ", fname.c_str()); logFlush();
		matrix epsInvEvecs; diagMatrix epsInvEigs;
		matrix Khalf = pow(dagger_symmetrize(K), 0.5);
		matrix epsInv = eye(nColumns) + Khalf * Xext * Khalf;
		epsInv.diagonalize(epsInvEvecs, epsInvEigs); //epsInv (symmetrized)
		FILE* fp=fopen(fname.c_str(), "w");
		epsInvEigs.print(fp, "%.15f\n");
		fclose(fp);
		logPrintf("Done.\n"); logFlush();
		//Print head:
		size_t iGzero = 0;
		for(const vector3<int>& iG: basis.iGarr)
		{	if(!iG.length_squared()) break;
			iGzero++;
		}
		assert(iGzero < basis.nbasis);
		matrix V0(1, nColumns);
		for(int j=0; j<nColumns; j++)
			V0.set(0,j, V.data()[V.index(j,iGzero)]);
		V0 *= 1./sqrt(trace(V0*dagger(V0)).abs()); //normalize
		double epsInvHead = trace(V0 * epsInv * dagger(V0)).abs();
		logPrintf("\thead(epsInv): %lg   1/head(epsInv): %lg\n", epsInvHead, 1./epsInvHead);
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
	//G-vectors:
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(e.dump.getFilename("pol_Gvectors").c_str(), "w");
		for(const vector3<int>& iG: basis.iGarr)
			fprintf(fp, "%d %d %d\n", iG[0], iG[1], iG[2]);
		fclose(fp);
	}
	logPrintf("Done.\n");
	logFlush();
}

string Polarizability::dkFilename(int ik, string varName) const
{
	ostringstream ikOss; ikOss << ik;
	//Create a map of substitutions:
	std::map<string,string> subMap;
	subMap["$VAR"] = varName;
	subMap["$q"] = ikOss.str();
	//Apply the substitutions:
	string fname = dkFilenamePattern;
	for(auto sub: subMap)
	{	size_t pos = fname.find(sub.first);
		if(pos != string::npos)
			fname.replace(pos, sub.first.length(), sub.second);
	}
	return fname;
}
