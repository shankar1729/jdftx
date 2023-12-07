/*-------------------------------------------------------------------
Copyright 2022 Brandon Li

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

#include "PerturbationInfo.h"
#include <core/Random.h>
#include <core/ScalarFieldIO.h>
#include <electronic/Everything.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecInfo.h>
#include <electronic/ElecMinimizer.h>
#include <perturb/PerturbationSolver.h>
#include <core/LatticeUtils.h>


Perturbation::Perturbation(const Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo) {};

AtomPerturbation::AtomPerturbation (unsigned int sp, unsigned int at, int iDir, const Everything& e) : Perturbation(e)
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

AtomPerturbation::AtomPerturbation (unsigned int sp, unsigned int at, vector3<> dirCartesian, const Everything& e) : Perturbation(e)
{
	AtomicMode m;
	
	m.sp = sp;
	m.at = at;
	m.dirCartesian = dirCartesian;

	m.dirLattice = e.gInfo.invR*m.dirCartesian;
	mode = m;
	init(e, e.eVars, e.eInfo);
}

void AtomPerturbation::init(const Everything &e, const ElecVars& eVars, const ElecInfo& eInfo)
{
	   Vatom.resize(eInfo.nStates);
	   dVatom.resize(eInfo.nStates);
	   VdagCatom.resize(eInfo.nStates);
	   dVdagCatom.resize(eInfo.nStates);
	::init(dCatom, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	dnatom.resize(eVars.Vscloc.size());
	nullToZero(dnatom, e.gInfo);
}

bool AtomPerturbation::isUltrasoft(const IonInfo& iInfo) { return iInfo.species[mode.sp]->isUltrasoft(); } //!< Does this atom use an ultrasoft potential

void PerturbationInfo::setup(const Everything &e, const ElecVars &eVars)
{
	solverParams.fpLog = globalLog;
	solverParams.linePrefix = "PerturbMinimize: ";
	if (!solverParams.nIterations)
		return;

	const ElecInfo &eInfo = e.eInfo;

	if (!commensurate)
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
		CdagdHC.resize(eInfo.nStates);
		CdagdHCatom.resize(eInfo.nStates);
		
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
		dnpq.resize(eVars.n.size());
		dnmq.resize(eVars.n.size());
		nullToZero(dnpq, e.gInfo);
		nullToZero(dnmq, e.gInfo);
	}
	
	if(dVext && dVext->dVexternalFilename.size())
	{
		if(commensurate) {
			dVext->dVext.resize(1);
			dVext->dVext[0] = ScalarFieldData::alloc(e.gInfo);
			logPrintf("Reading external perturbation potential from '%s'\n", dVext->dVexternalFilename[0].c_str());
			loadRawBinary(dVext->dVext[0], dVext->dVexternalFilename[0].c_str());
		} else {
			dVext->dVextpq.resize(1);
			dVext->dVextmq.resize(1);
			dVext->dVextpq[0] = complexScalarFieldData::alloc(e.gInfo);
			logPrintf("Reading external perturbation potential from '%s'\n", dVext->dVexternalFilename[0].c_str());
			loadRawBinary(dVext->dVextpq[0], dVext->dVexternalFilename[0].c_str());
			
			dVext->dVextmq = conj(dVext->dVextpq);
		}
	}


	if(drhoExt && drhoExt->drhoExtFilename.size())
	{
		if(commensurate) {
			ScalarField temp(ScalarFieldData::alloc(e.gInfo));
			logPrintf("Reading external charge perturbation from '%s'\n", drhoExt->drhoExtFilename[0].c_str());
			loadRawBinary(temp, drhoExt->drhoExtFilename[0].c_str());
			drhoExt->drhoExt = J(temp);
		} else {
			complexScalarField temp(complexScalarFieldData::alloc(e.gInfo));
			logPrintf("Reading external charge perturbation from '%s'\n", drhoExt->drhoExtFilename[0].c_str());
			loadRawBinary(temp, drhoExt->drhoExtFilename[0].c_str());

			drhoExt->drhoExtpq = J(temp);
			drhoExt->drhoExtmq = J(conj(I(drhoExt->drhoExtpq)));
		}
	}

	if (!commensurate)
		read(e.eInfo, e, Cinc, wfnsFilename.c_str(), nullptr);
	
	checkSupportedFeatures(e, e.eInfo);
}


void PerturbationInfo::read(const ElecInfo &eInfo, const Everything &e, std::vector<ColumnBundle>& C, const char *fname, const ElecInfo::ColumnBundleReadConversion* conversion) const
{	if(conversion && conversion->realSpace)
	{
		die("Real space incommensurate wavefunctions are required.\n")
	}
	else
	{	//Check if a conversion is actually needed:
		std::vector<ColumnBundle> Ctmp_kplusq(eInfo.qStop);
		std::vector<ColumnBundle> Ctmp_kminusq(eInfo.qStop);
		std::vector<Basis> basisTmp_kplusq(eInfo.qStop);
		std::vector<Basis> basisTmp_kminusq(eInfo.qStop);
		std::vector<long> nBytes_kplusq(mpiWorld->nProcesses(), 0); //total bytes to be read on each process
		std::vector<long> nBytes_kminusq(mpiWorld->nProcesses(), 0);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int kpq = q;
			int kmq = q + eInfo.nStates;
			bool needTmp = false, customBasis = false;
			int nCols = C[kpq].nCols();
			if(conversion)
			{
				double EcutOld = conversion->EcutOld ? conversion->EcutOld : conversion->Ecut;
				customBasis = (EcutOld!=conversion->Ecut);
				if(customBasis)
				{	needTmp = true;
					logSuspend();
					basisTmp_kplusq[q].setup(*(C[kpq].basis->gInfo), *(C[kpq].basis->iInfo), EcutOld, C[kpq].qnum->k);
					basisTmp_kminusq[q].setup(*(C[kmq].basis->gInfo), *(C[kmq].basis->iInfo), EcutOld, C[kmq].qnum->k);
					logResume();
				}
			}
			const Basis* basis = customBasis ? &basisTmp_kplusq[q] : C[kpq].basis;
			int nSpinor = C[kpq].spinorLength();
			if(needTmp) Ctmp_kplusq[q].init(nCols, basis->nbasis*nSpinor, basis, C[kpq].qnum);
			nBytes_kplusq[mpiWorld->iProcess()] += nCols * basis->nbasis*nSpinor * sizeof(complex);

			basis = customBasis ? &basisTmp_kminusq[q] : C[kmq].basis;
			nSpinor = C[kmq].spinorLength();
			if(needTmp) Ctmp_kminusq[q].init(nCols, basis->nbasis*nSpinor, basis, C[kmq].qnum);
			nBytes_kminusq[mpiWorld->iProcess()] += nCols * basis->nbasis*nSpinor * sizeof(complex);
		}
		//Sync nBytes:
		if(mpiWorld->nProcesses()>1)
			for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++) {
				mpiWorld->bcast( nBytes_kplusq[iSrc], iSrc);
				mpiWorld->bcast( nBytes_kminusq[iSrc], iSrc);
			}
		//Compute offset of current process, and expected file length:
		long offset_kpq=0, fsize=0;
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
		{	if(iSrc<mpiWorld->iProcess()) offset_kpq += nBytes_kplusq[iSrc];
			fsize += nBytes_kplusq[iSrc];
		}

		long offset_kmq=fsize;
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
		{	if(iSrc<mpiWorld->iProcess()) offset_kmq += nBytes_kminusq[iSrc];
			fsize += nBytes_kminusq[iSrc];
		}

		//Read data into Ytmp or Y as appropriate, and convert if necessary:
		MPIUtil::File fp; mpiWorld->fopenRead(fp, fname, fsize, "Hint: Did you specify the correct EcutOld?\n");
		mpiWorld->fseek(fp, offset_kpq, SEEK_SET);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int kpq = q;
			ColumnBundle& Ycur = Ctmp_kplusq[q] ? Ctmp_kplusq[q] : C[kpq];
			mpiWorld->freadData(Ycur, fp);
			if(Ctmp_kplusq[q]) //apply conversions:
			{	if(Ctmp_kplusq[q].basis!=C[kpq].basis)
				{	int nSpinor = C[kpq].spinorLength();
					for(int b=0; b<std::min(C[kpq].nCols(), Ctmp_kplusq[q].nCols()); b++)
						for(int s=0; s<nSpinor; s++)
							C[kpq].setColumn(b,s, Ctmp_kplusq[q].getColumn(b,s)); //convert using the full G-space as an intermediate
				}
				else
				{	if(Ctmp_kplusq[q].nCols()<C[kpq].nCols()) C[kpq].setSub(0, Ctmp_kplusq[q]);
					else C[kpq] = Ctmp_kplusq[q].getSub(0, C[kpq].nCols());
				}
				Ctmp_kplusq[q].free();
			}
		}

		mpiWorld->fseek(fp, offset_kmq, SEEK_SET);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			int kmq = q + eInfo.nStates;
			ColumnBundle& Ycur = Ctmp_kminusq[q] ? Ctmp_kminusq[q] : C[kmq];
			mpiWorld->freadData(Ycur, fp);
			if(Ctmp_kminusq[q]) //apply conversions:
			{	if(Ctmp_kminusq[q].basis!=C[kmq].basis)
				{	int nSpinor = C[kmq].spinorLength();
					for(int b=0; b<std::min(C[kmq].nCols(), Ctmp_kminusq[q].nCols()); b++)
						for(int s=0; s<nSpinor; s++)
							C[kmq].setColumn(b,s, Ctmp_kminusq[q].getColumn(b,s)); //convert using the full G-space as an intermediate
				}
				else
				{	if(Ctmp_kminusq[q].nCols()<C[kmq].nCols()) C[kmq].setSub(0, Ctmp_kminusq[q]);
					else C[kmq] = Ctmp_kminusq[q].getSub(0, C[kmq].nCols());
				}
				Ctmp_kminusq[q].free();
			}
		}
		mpiWorld->fclose(fp);

		logPrintf("Successfully read band minimized wavefunctions.\n");
	}
}

void PerturbationInfo::setupkpoints(const Everything &e, const ElecInfo &eInfo)
{
	kplusq_vectors = eInfo.qnums;
	kminusq_vectors = eInfo.qnums;

	for(int q=0; q<eInfo.nStates; q++) {
		kplusq_vectors[q].k = kplusq_vectors[q].k + qvec;
		kminusq_vectors[q].k = kminusq_vectors[q].k - qvec;
	}

	const GridInfo& gInfoBasis = e.gInfoWfns ? *e.gInfoWfns : e.gInfo;

	if(e.cntrl.basisKdep != BasisKpointDep)
		die("Please use k-point dependent basis.");

	kplusq_basis.resize(eInfo.nStates);
	kminusq_basis.resize(eInfo.nStates);

	for(int q=0; q<eInfo.nStates; q++)
	{
		kplusq_basis[q].setup(gInfoBasis, e.iInfo, e.cntrl.Ecut, kplusq_vectors[q].k);
		kminusq_basis[q].setup(gInfoBasis, e.iInfo, e.cntrl.Ecut, kminusq_vectors[q].k);
	}

	logPrintf("Printing k, k+q, and k-q vectors:\n");
	for(int q=0; q<eInfo.nStates; q++)
	{
		logPrintf("index: %5d, k=[ %+.7f %+.7f %+.7f ]  weight=%.9f\n", q, eInfo.qnums[q].k[0], eInfo.qnums[q].k[1], eInfo.qnums[q].k[2], eInfo.qnums[q].weight);
		logPrintf("index: %5d, k+q=[ %+.7f %+.7f %+.7f ]  weight=%.9f\n", q, kplusq_vectors[q].k[0], kplusq_vectors[q].k[1], kplusq_vectors[q].k[2], kplusq_vectors[q].weight);
		logPrintf("index: %5d, k-q=[ %+.7f %+.7f %+.7f ]  weight=%.9f\n", q, kminusq_vectors[q].k[0], kminusq_vectors[q].k[1], kminusq_vectors[q].k[2], kminusq_vectors[q].weight);
	}
}

void PerturbationInfo::initInc(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const ElecInfo* eInfo)
{
	Y.resize(nbundles);
	if(ncols && eInfo)
	{	assert(nbundles == 2*eInfo->qStop);
		for(int q=eInfo->qStart; q<eInfo->qStop; q++) {
			int kpq = q ;
			int kmq = q + eInfo->nStates;
			Y[kpq].init(ncols, kplusq_basis[q].nbasis * eInfo->spinorLength(), &kplusq_basis[q], &kplusq_vectors[q], isGpuEnabled());
			Y[kmq].init(ncols, kminusq_basis[q].nbasis * eInfo->spinorLength(), &kminusq_basis[q], &kminusq_vectors[q], isGpuEnabled());
		}
	}
}

void PerturbationInfo::checkSupportedFeatures(const Everything &e, const ElecInfo &eInfo)
{
	if(!(eInfo.fillingsUpdate==ElecInfo::FillingsConst && eInfo.scalarFillings))
		die("Constant fillings only.\n");
	if(e.exCorr.exxFactor())
		die("Variational perturbation currently does not support exact exchange.\n");
	if(e.eInfo.hasU)
		die("Variational perturbation currently does not support DFT+U.\n");
	if (e.exCorr.needsKEdensity())
		die("Variational perturbation currently does not support KE density dependent functionals.\n");
	if(!e.exCorr.hasEnergy())
		die("Variational perturbation does not support potential functionals.\n");
	if(e.exCorr.orbitalDep)
		die("Variational perturbation currently does not support orbital dependent potential functionals.\n");
	if(e.eVars.fluidParams.fluidType != FluidNone)
		die("Variational perturbation does not support fluids.\n");
	if (e.eVars.n.size() > 1)
		die("Multiple spin channels not supported yet.\n");
	if (e.coulombParams.geometry != CoulombParams::Geometry::Periodic && !commensurate)
		die("Periodic coloumb interaction required for incommensurate perturbations.\n")
		
	if (e.symm.mode != SymmetriesNone)
	{
		if (!commensurate)
			die("Symmetries are not supported with incommensurate perturbations.\n");
		
		if (dVext) {
			ScalarFieldArray dVext_sym = clone(dVext->dVext);
			e.symm.symmetrize(dVext_sym);
			if (nrm2(dVext_sym[0] - dVext->dVext[0])/nrm2(dVext->dVext[0]) < symmThreshold)
			logPrintf("Warning: dVext->dVext does not obey symmetries.\n");
		}
	}
		
	for (auto sp: e.iInfo.species) {
		if (sp->isRelativistic())
			die("Relativistic potentials are not supported yet.\n");
		if (sp->isUltrasoft() & !commensurate)
			die("Ultrasoft potentials are compatible with commensurate perturbations only.\n");
	}

	if (e.exCorr.needFiniteDifferencing())
	logPrintf("Excorr analytical derivative not available. Using finite differencing instead.\n");

	if (!commensurate) {
		if(e.eVars.wfnsFilename.empty() && wfnsFilename.empty()) {
			die("Currently, incommensurate perturbations require ground state wavefunctions to be loaded in.\n");
			//computeIncommensurateWfns();
		} else if (!e.eVars.wfnsFilename.empty() && wfnsFilename.empty()) {
			die("Please specify offset wavefunctions.\n")
		} else if (e.eVars.wfnsFilename.empty() && !wfnsFilename.empty()) {
			die("Please specify ground state wavefunctions.\n")
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


void PerturbationInfo::sampleCB (ColumnBundle C, std::string name) {
	if (C.nData() < 10) return;
	double *dat = (double*)(C.getColumn(0, 0)->dataPref());
	logPrintf("ColumnBundle %s values: %g %g %g %g %g %g %g %g %g %g\n", name.c_str(), dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void PerturbationInfo::sampleMat (matrix C, std::string name) {
	if (C.nRows() < 3 || C.nCols() < 3) return;
	logPrintf("Matrix %s values %g %g %g %g %g %g %g %g %g\n", name.c_str(), C(0,0).x, C(0,1).x, C(0,2).x, C(1,0).x,C(1,1).x,C(1,2).x,C(2,0).x,C(2,1).x,C(2,2).x);
}

void PerturbationInfo::sampleField (ScalarField V, std::string name) {
	if (V->nElem < 10) return;
	double *dat = V->dataPref();
	logPrintf("ScalarField %s values %g %g %g %g %g %g %g %g %g %g\n", name.c_str(), dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void PerturbationInfo::sampleField (ScalarFieldTilde V, std::string name) {
	if (V->nElem < 10) return;
	double *dat = (double*)(V->dataPref());
	logPrintf("ScalarFieldTidle %s values %g %g %g %g %g %g %g %g %g %g\n", name.c_str(), dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}
