/*
 * PerturbationInfo.cpp
 *
 *  Created on: Jul 25, 2022
 *      Author: brandon
 */

#include "PerturbationInfo.h"
#include <core/Random.h>
#include <core/ScalarFieldIO.h>
#include <electronic/Everything.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecInfo.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/PerturbationSolver.h>

PerturbationInfo::PerturbationInfo() {
}

PerturbationInfo::~PerturbationInfo() {
	// TODO Auto-generated destructor stub
}

void PerturbationInfo::setup(const Everything &e, const ElecVars &eVars) {

	if (!e.vptParams.nIterations)
		return;

	const ElecInfo &eInfo = e.eInfo;

	setupkpoints(e, eInfo);

	//Allocate variables

	if (!incommensurate) {
		init(dY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(HdC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dHC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(HdCatom, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dHCatom, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dGradPsi, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dGradTau, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(HC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(dCatom, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(OC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		dU.resize(eInfo.nStates);
		dUsqrtinvatom.resize(eInfo.nStates);
		dVscloc.resize(eVars.Vscloc.size());
		dVsclocatom.resize(eVars.Vscloc.size());
		dnatom.resize(eVars.Vscloc.size());
		Vatom_cached.resize(eInfo.nStates);
		dVatom_cached.resize(eInfo.nStates);
		VdagCatom_cached.resize(eInfo.nStates);
		dVdagCatom_cached.resize(eInfo.nStates);
	} else {
		initInc(dY, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		initInc(dC, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		initInc(dGradPsi, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		initInc(dGradTau, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		initInc(Cinc, 2*eInfo.nStates, eInfo.nBands, &eInfo);
        dVsclocpq.resize(eVars.Vscloc.size());
        dVsclocmq.resize(eVars.Vscloc.size());
	}
	
	if(dVexternalFilename.size())
	{
		if(!incommensurate)
		{
			dVext.resize(dVexternalFilename.size());
			for(unsigned s=0; s<dVext.size(); s++)
			{	dVext[s] = ScalarFieldData::alloc(e.gInfo);
				logPrintf("Reading external perturbation potential from '%s'\n", dVexternalFilename[s].c_str());
				loadRawBinary(dVext[s], dVexternalFilename[s].c_str());
			}
			if(dVext.size()==1 && eVars.n.size()==2) //Replicate potential for second spin:
				dVext.push_back(dVext[0]->clone());
			
		} else {
			dVextpq.resize(dVexternalFilename.size());
			for(unsigned s=0; s<dVextpq.size(); s++)
			{	dVextpq[s] = complexScalarFieldData::alloc(e.gInfo);
				logPrintf("Reading external perturbation potential from '%s'\n", dVexternalFilename[s].c_str());
				loadRawBinary(dVextpq[s], dVexternalFilename[s].c_str());
			}
			if(dVextpq.size()==1 && eVars.n.size()==2) //Replicate potential for second spin:
				dVextpq.push_back(dVextpq[0]->clone());
			
			dVextmq = conj(dVextpq);
		}
        VextPerturbationExists = true;
	} else {
		if(!incommensurate) {
			dVext.resize(eVars.Vexternal.size());
			nullToZero(dVext, e.gInfo);
		} else {
			dVextpq.resize(eVars.Vexternal.size());
			dVextmq.resize(eVars.Vexternal.size());
			nullToZero(dVextpq, e.gInfo);
			nullToZero(dVextmq, e.gInfo);
		}
        VextPerturbationExists = false;
		logPrintf("Note: dVext not loaded. Setting external perturbation to zero.\n");
	}
	//assert(dVext.size() == eVars.Vexternal.size());

	if (incommensurate)
		read(e.eInfo, e, Cinc, wfnsFilename.c_str(), nullptr);

	dn.resize(eVars.n.size());
	dnpq.resize(eVars.n.size());
	dnmq.resize(eVars.n.size());
}


void PerturbationInfo::read(const ElecInfo &eInfo, const Everything &e, std::vector<ColumnBundle>& C, const char *fname, const ElecInfo::ColumnBundleReadConversion* conversion) const
{	if(conversion && conversion->realSpace)
	{
		die("Please use fourier space wavefunctions")
	}
	else
	{	//Check if a conversion is actually needed:
		std::vector<ColumnBundle> Ytmp_Tk(eInfo.qStop);
		std::vector<ColumnBundle> Ytmp_Tinvk(eInfo.qStop);
		std::vector<Basis> basisTmp_Tk(eInfo.qStop);
		std::vector<Basis> basisTmp_Tinvk(eInfo.qStop);
		std::vector<long> nBytes_Tk(mpiWorld->nProcesses(), 0); //total bytes to be read on each process
		std::vector<long> nBytes_Tinvk(mpiWorld->nProcesses(), 0);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int Tk = q ;
			int Tinvk = q + eInfo.nStates;
			bool needTmp = false, customBasis = false;
			int nCols = C[Tk].nCols();
			if(conversion)
			{	if(conversion->nBandsOld && conversion->nBandsOld!=nCols)
				{	nCols = conversion->nBandsOld;
					needTmp = true;
				}
				double EcutOld = conversion->EcutOld ? conversion->EcutOld : conversion->Ecut;
				customBasis = (EcutOld!=conversion->Ecut);
				if(customBasis)
				{	needTmp = true;
					logSuspend();
					basisTmp_Tk[q].setup(*(C[Tk].basis->gInfo), *(C[Tk].basis->iInfo), EcutOld, C[Tk].qnum->k);
					basisTmp_Tinvk[q].setup(*(C[Tinvk].basis->gInfo), *(C[Tinvk].basis->iInfo), EcutOld, C[Tinvk].qnum->k);
					logResume();
				}
			}
			const Basis* basis = customBasis ? &basisTmp_Tk[q] : C[Tk].basis;
			int nSpinor = C[Tk].spinorLength();
			if(needTmp) Ytmp_Tk[q].init(nCols, basis->nbasis*nSpinor, basis, C[Tk].qnum);
			nBytes_Tk[mpiWorld->iProcess()] += nCols * basis->nbasis*nSpinor * sizeof(complex);

			basis = customBasis ? &basisTmp_Tinvk[q] : C[Tinvk].basis;
			nSpinor = C[Tinvk].spinorLength();
			if(needTmp) Ytmp_Tinvk[q].init(nCols, basis->nbasis*nSpinor, basis, C[Tinvk].qnum);
			nBytes_Tinvk[mpiWorld->iProcess()] += nCols * basis->nbasis*nSpinor * sizeof(complex);
		}
		//Sync nBytes:
		if(mpiWorld->nProcesses()>1)
			for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++) {
				mpiWorld->bcast(nBytes_Tk[iSrc], iSrc);
				mpiWorld->bcast(nBytes_Tinvk[iSrc], iSrc);
			}
		//Compute offset of current process, and expected file length:
		long offset_Tk=0, fsize=0;
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
		{	if(iSrc<mpiWorld->iProcess()) offset_Tk += nBytes_Tk[iSrc];
			fsize += nBytes_Tk[iSrc];
		}

		long offset_Tinvk=fsize;
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
		{	if(iSrc<mpiWorld->iProcess()) offset_Tinvk += nBytes_Tinvk[iSrc];
			fsize += nBytes_Tinvk[iSrc];
		}

		//Read data into Ytmp or Y as appropriate, and convert if necessary:
		MPIUtil::File fp; mpiWorld->fopenRead(fp, fname, fsize,
			(e.vibrations and eInfo.qnums.size()>1)
			? "Hint: Vibrations requires wavefunctions without symmetries:\n"
				"either don't read in state, or consider using phonon instead.\n"
			: "Hint: Did you specify the correct nBandsOld, EcutOld and kdepOld?\n");
		mpiWorld->fseek(fp, offset_Tk, SEEK_SET);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int Tk = q ;
			ColumnBundle& Ycur = Ytmp_Tk[q] ? Ytmp_Tk[q] : C[Tk];
			mpiWorld->freadData(Ycur, fp);
			if(Ytmp_Tk[q]) //apply conversions:
			{	if(Ytmp_Tk[q].basis!=C[Tk].basis)
				{	int nSpinor = C[Tk].spinorLength();
					for(int b=0; b<std::min(C[Tk].nCols(), Ytmp_Tk[q].nCols()); b++)
						for(int s=0; s<nSpinor; s++)
							C[Tk].setColumn(b,s, Ytmp_Tk[q].getColumn(b,s)); //convert using the full G-space as an intermediate
				}
				else
				{	if(Ytmp_Tk[q].nCols()<C[Tk].nCols()) C[Tk].setSub(0, Ytmp_Tk[q]);
					else C[Tk] = Ytmp_Tk[q].getSub(0, C[Tk].nCols());
				}
				Ytmp_Tk[q].free();
			}
		}

		mpiWorld->fseek(fp, offset_Tinvk, SEEK_SET);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int Tinvk = q + eInfo.nStates;
			ColumnBundle& Ycur = Ytmp_Tinvk[q] ? Ytmp_Tinvk[q] : C[Tinvk];
			mpiWorld->freadData(Ycur, fp);
			if(Ytmp_Tinvk[q]) //apply conversions:
			{	if(Ytmp_Tinvk[q].basis!=C[Tinvk].basis)
				{	int nSpinor = C[Tinvk].spinorLength();
					for(int b=0; b<std::min(C[Tinvk].nCols(), Ytmp_Tinvk[q].nCols()); b++)
						for(int s=0; s<nSpinor; s++)
							C[Tinvk].setColumn(b,s, Ytmp_Tinvk[q].getColumn(b,s)); //convert using the full G-space as an intermediate
				}
				else
				{	if(Ytmp_Tinvk[q].nCols()<C[Tinvk].nCols()) C[Tinvk].setSub(0, Ytmp_Tinvk[q]);
					else C[Tinvk] = Ytmp_Tinvk[q].getSub(0, C[Tinvk].nCols());
				}
				Ytmp_Tinvk[q].free();
			}
		}
		mpiWorld->fclose(fp);

		logPrintf("Successfully read band minimized wavefunctions");
	}
}

void PerturbationInfo::setupkpoints(const Everything &e, const ElecInfo &eInfo)
{
	Tk_vectors = eInfo.qnums;
	Tinvk_vectors = eInfo.qnums;

	for(int q=0; q<eInfo.nStates; q++) {
		Tk_vectors[q].k = Tk_vectors[q].k + qvec;
		Tinvk_vectors[q].k = Tinvk_vectors[q].k - qvec;
	}


	const GridInfo& gInfoBasis = e.gInfoWfns ? *e.gInfoWfns : e.gInfo;

	if(e.cntrl.basisKdep != BasisKpointDep)
		die("Please use k-point dependent basis.");

	Tk_basis.resize(eInfo.nStates);
	Tinvk_basis.resize(eInfo.nStates);

	for(int q=0; q<eInfo.nStates; q++)
	{
		//TODO check Ecut is same
		Tk_basis[q].setup(gInfoBasis, e.iInfo, e.cntrl.Ecut, Tk_vectors[q].k);
		Tinvk_basis[q].setup(gInfoBasis, e.iInfo, e.cntrl.Ecut, Tinvk_vectors[q].k);
	}

	logPrintf("Printing k, k+q , and k-q vectors\n");
	for(int q=0; q<eInfo.nStates; q++)
	{
		QuantumNumber qnum = eInfo.qnums[q];
		logPrintf("%5d  [ %+.7f %+.7f %+.7f ]  %.9f\n", q, qnum.k[0], qnum.k[1], qnum.k[2], qnum.weight);
		qnum = Tk_vectors[q];
		logPrintf("%5d  [ %+.7f %+.7f %+.7f ]  %.9f\n", q, qnum.k[0], qnum.k[1], qnum.k[2], qnum.weight);
		qnum = Tinvk_vectors[q];
		logPrintf("%5d  [ %+.7f %+.7f %+.7f ]  %.9f\n", q, qnum.k[0], qnum.k[1], qnum.k[2], qnum.weight);
	}
}

void PerturbationInfo::initInc(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const ElecInfo* eInfo)
{
	Y.resize(nbundles);
	if(ncols && eInfo)
	{	assert(nbundles == 2*eInfo->qStop);
		for(int q=eInfo->qStart; q<eInfo->qStop; q++) {
			int Tk = q ;
			int Tinvk = q + eInfo->nStates;
			Y[Tk].init(ncols, Tk_basis[q].nbasis * eInfo->spinorLength(), &Tk_basis[q], &Tk_vectors[q], isGpuEnabled());
			Y[Tinvk].init(ncols, Tinvk_basis[q].nbasis * eInfo->spinorLength(), &Tinvk_basis[q], &Tinvk_vectors[q], isGpuEnabled());
		}
	}
}

void PerturbationInfo::updateExcorrCache(const ExCorr& exc, const GridInfo& gInfo, const ScalarField& n)
{
	int iDirStart, iDirStop;
	TaskDivision(3, mpiWorld).myRange(iDirStart, iDirStop);
	{
		for(int i=iDirStart; i<iDirStop; i++) {
			IDJn_cached[i] = I(D(J(n),i));
        }

		nullToZero(sigma_cached, gInfo);
		initZero(sigma_cached);

		for(int i=iDirStart; i<iDirStop; i++)
			sigma_cached += IDJn_cached[i] * IDJn_cached[i];
		sigma_cached->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	
	exc.getSecondDerivatives(n, e_nn_cached, e_sigma_cached, e_nsigma_cached, e_sigmasigma_cached, 1e-9, &sigma_cached);
}
