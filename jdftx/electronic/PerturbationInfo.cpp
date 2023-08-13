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
#include <core/LatticeUtils.h>


AtomPerturbation::AtomPerturbation (int sp, int at, int iDir, const Everything& e)
{
    AtomicMode m;
    vector3<> dirCartesian;

    dirCartesian[iDir] = 1;
    m.sp = sp;
    m.at = at;
    m.dirCartesian = dirCartesian;

    m.dirLattice = e.gInfo.invR*m.dirCartesian;
    mode = m;
    init(e, e.eVars, e.eInfo);
}

AtomPerturbation::AtomPerturbation (int sp, int at, vector3<> dirCartesian, const Everything& e)
{
    AtomicMode m;
	
    m.sp = sp;
    m.at = at;
    m.dirCartesian = dirCartesian;

    m.dirLattice = e.gInfo.invR*m.dirCartesian;
    mode = m;
    init(e, e.eVars, e.eInfo);
}

bool AtomPerturbation::sameAtom(const std::shared_ptr<AtomPerturbation> pert) {
	return pert->mode.sp == mode.sp && pert->mode.at == mode.at;
}

void AtomPerturbation::init(const Everything &e, const ElecVars& eVars, const ElecInfo& eInfo)
{
	Vatom_cached.resize(eInfo.nStates);
	dVatom_cached.resize(eInfo.nStates);
	VdagCatom_cached.resize(eInfo.nStates);
	dVdagCatom_cached.resize(eInfo.nStates);
	::init(dCatom, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	dnatom.resize(eVars.Vscloc.size());
	nullToZero(dnatom, e.gInfo);
}

bool AtomPerturbation::isUltrasoft(const IonInfo& iInfo) {
	return iInfo.species[mode.sp]->isUltrasoft();
}

PerturbationInfo::PerturbationInfo() : dVext(0), datom(0) {
	
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

	if (commensurate) {
		init(dC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(HC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(OC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		init(grad, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		dU.resize(eInfo.nStates);
		dUmhalfatom.resize(eInfo.nStates);
		dHsub.resize(eInfo.nStates);
		dHsubatom.resize(eInfo.nStates);
		dVscloc.resize(eVars.Vscloc.size());
		dVsclocTau.resize(eVars.Vscloc.size());
		dn.resize(eVars.Vscloc.size());
		nullToZero(dn, e.gInfo);
		if (datom) datom->init(e, eVars, eInfo);
	} else {
		initInc(dC, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		initInc(Cinc, 2*eInfo.nStates, eInfo.nBands, &eInfo);
        dVsclocpq.resize(eVars.Vscloc.size());
        dVsclocmq.resize(eVars.Vscloc.size());
	}
	
	if(dVext && dVext->dVexternalFilename.size())
	{
		if(commensurate)
		{
			dVext->dVext.resize(dVext->dVexternalFilename.size());
			for(unsigned s=0; s<dVext->dVext.size(); s++)
			{	dVext->dVext[s] = ScalarFieldData::alloc(e.gInfo);
				logPrintf("Reading external perturbation potential from '%s'\n", dVext->dVexternalFilename[s].c_str());
				loadRawBinary(dVext->dVext[s], dVext->dVexternalFilename[s].c_str());
			}
			if(dVext->dVext.size()==1 && eVars.n.size()==2) //Replicate potential for second spin:
				dVext->dVext.push_back(dVext->dVext[0]->clone());
			
		} else {
			dVext->dVextpq.resize(dVext->dVexternalFilename.size());
			for(unsigned s=0; s<dVext->dVextpq.size(); s++)
			{	dVext->dVextpq[s] = complexScalarFieldData::alloc(e.gInfo);
				logPrintf("Reading external perturbation potential from '%s'\n", dVext->dVexternalFilename[s].c_str());
				loadRawBinary(dVext->dVextpq[s], dVext->dVexternalFilename[s].c_str());
			}
			if(dVext->dVextpq.size()==1 && eVars.n.size()==2) //Replicate potential for second spin:
				dVext->dVextpq.push_back(dVext->dVextpq[0]->clone());
			
			dVext->dVextmq = conj(dVext->dVextpq);
		}
	}
	//assert(dVext->dVext.size() == eVars.Vexternal.size());

	if (!commensurate)
		read(e.eInfo, e, Cinc, wfnsFilename.c_str(), nullptr);

	dn.resize(eVars.n.size());
	dnpq.resize(eVars.n.size());
	dnmq.resize(eVars.n.size());
	
	checkSupportedFeatures(e, e.eInfo);
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

void PerturbationInfo::checkSupportedFeatures(const Everything &e, const ElecInfo &eInfo)
{
	if(!(eInfo.fillingsUpdate==ElecInfo::FillingsConst && eInfo.scalarFillings))
		die("Constant fillings only.");
	if(e.exCorr.exxFactor())
		die("Variational perturbation currently does not support exact exchange.");
	if(e.eInfo.hasU)
		die("Variational perturbation currently does not support DFT+U.");
	if (e.exCorr.needsKEdensity())
		die("Variational perturbation currently does not support KE density dependent functionals.");
	if(!e.exCorr.hasEnergy())
		die("Variational perturbation does not support potential functionals.");
	if(e.exCorr.orbitalDep)
		die("Variational perturbation currently does not support orbital dependent potential functionals.");
	if(e.eVars.fluidParams.fluidType != FluidNone)
		die("Variational perturbation does not support fluids.");
	if (e.eVars.n.size() > 1)
		die("Multiple spin channels not supported yet.");
	if (e.coulombParams.geometry != CoulombParams::Geometry::Periodic && !commensurate)
		die("Periodic coloumb interaction required for incommensurate perturbations.")

		
		
	if (e.symm.mode != SymmetriesNone)
	{
		if (!commensurate)
			die("Symmetries are not supported with incommensurate perturbations.");
		
		logPrintf("Warning: VPT has not been tested with symmetries");
		if (dVext) {
			ScalarFieldArray dVext_sym = clone(dVext->dVext);
			e.symm.symmetrize(dVext_sym);
			if (nrm2(dVext_sym[0] - dVext->dVext[0])/nrm2(dVext->dVext[0]) < symmThreshold)
			logPrintf("Warning: dVext->dVext does not obey symmetries.");
		}
	}
		
	for (auto sp: e.iInfo.species) {
		if (sp->isRelativistic())
			die("Relativistic potentials are not supported yet.");
		if (sp->isUltrasoft() & !commensurate)
			die("Ultrasoft potentials are compatible with commensurate perturbations only.");
	}

	if (e.exCorr.needFiniteDifferencing())
	logPrintf("Excorr analytical derivative not available. Using finite differencing instead.\n");

	if (!commensurate) {
		if(e.eVars.wfnsFilename.empty() && wfnsFilename.empty()) {
			die("Currently, incommensurate perturbations require ground state wavefunctions to be loaded in.\n");
			//computeIncommensurateWfns();
		} else if (!e.eVars.wfnsFilename.empty() && wfnsFilename.empty()) {
			die("Please specify offset wavefunctions")
		} else if (e.eVars.wfnsFilename.empty() && !wfnsFilename.empty()) {
			die("Please specify ground state wavefunctions")
		}
	}
}

bool PerturbationInfo::densityAugRequired(const Everything &e) {
	for (auto sp: e.iInfo.species) {
		if (sp->isUltrasoft())
			return true;
	}
	return false;
}
