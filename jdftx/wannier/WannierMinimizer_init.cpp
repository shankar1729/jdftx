/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <wannier/WannierMinimizer.h>
#include <core/WignerSeitz.h>
#include <core/ScalarFieldIO.h>
#include <core/Blip.h>

WannierMinimizer::WannierMinimizer(const Everything& e, const Wannier& wannier, bool needSuperOverride)
: e(e), wannier(wannier), sym(e.symm.getMatrices()),
	nCenters(wannier.nCenters), nFrozen(wannier.nFrozen), nBands(e.eInfo.nBands),
	nSpins(e.eInfo.nSpins()), qCount(e.eInfo.qnums.size()/nSpins),
	nSpinor(e.eInfo.spinorLength()),
	rSqExpect(nCenters), rExpect(nCenters), pinned(nCenters, false), rPinned(nCenters),
	needSuper(needSuperOverride || wannier.saveWfns || wannier.saveWfnsRealSpace || wannier.numericalOrbitalsFilename.length()),
	nPhononModes(0)
{
	//Create supercell grid:
	logPrintf("\n---------- Initializing supercell grid for Wannier functions ----------\n");
	const Supercell& supercell = *(e.coulombParams.supercell);
	gInfoSuper.R = supercell.Rsuper;
	gInfoSuper.GmaxRho = e.gInfo.Gmax;
	if(needSuper) gInfoSuper.initialize(true, sym);

	//Determine and output the band ranges:
	for(int iSpin: wannier.iSpinArr)
	{	std::vector<double> eMin(nBands, DBL_MAX), eMax(nBands, -DBL_MAX);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) if(e.eInfo.qnums[q].index()==iSpin)
			for(int b=0; b<nBands; b++)
			{	eMin[b] = std::min(eMin[b], e.eVars.Hsub_eigs[q][b]);
				eMax[b] = std::max(eMax[b], e.eVars.Hsub_eigs[q][b]);
			}
		mpiWorld->allReduceData(eMin, MPIUtil::ReduceMin);
		mpiWorld->allReduceData(eMax, MPIUtil::ReduceMax);
		if(mpiWorld->isHead())
		{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfBandRanges", &iSpin);
			logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
			FILE* fp = fopen(fname.c_str(), "w");
			for(int b=0; b<nBands; b++)
				fprintf(fp, "%+10.5lf %+10.5lf\n", eMin[b], eMax[b]);
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
	}
	
	//Create a list of kpoints:
	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	kMesh.resize(supercell.kmeshTransform.size());
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	KmeshEntry& ki = kMesh[ik];
		//Initialize kpoint:
		const Supercell::KmeshTransform& src = supercell.kmeshTransform[ik];
		//--- Copy over base class KMeshTransform:
		(Supercell::KmeshTransform&)ki.point = src;
		//--- Initialize base class QuantumNumber:
		ki.point.k = (~sym[src.iSym].rot) * (src.invert * qnums[src.iReduced].k) + src.offset;
		ki.point.weight = 1./kMesh.size();
		ki.point.spin = 0;
		kpoints.insert(ki.point);
	}
	
	//Determine overall Bloch wavevector of supercell (if any)
	if(needSuper)
	{	qnumSuper.weight = 1.;
		qnumSuper.spin = 0.;
		qnumSuper.k = kMesh[0].point.k * supercell.super;
		for(int l=0; l<3; l++)
			qnumSuper.k[l] -= floor(0.5+qnumSuper.k[l]);
		if(qnumSuper.k.length_squared()>symmThresholdSq)
		{	logPrintf("WARNING: k-mesh does not contain Gamma point. Orbitals will not be strictly periodic on supercell,\n"
				"\tbecause of an overall Bloch wave-vector: ");
			qnumSuper.k.print(globalLog, " %lf ");
		}
	}
	
	//Determine distribution amongst processes:
	kDivision.init(kMesh.size(), mpiWorld);
	kDivision.myRange(ikStart, ikStop);
	
	//Initialized pinned arrays:
	for(int n=nFrozen; n<nCenters; n++)
	{	pinned[n] = wannier.trialOrbitals[n-nFrozen].pinned;
		rPinned[n] = e.gInfo.R * wannier.trialOrbitals[n-nFrozen].xCenter;
	}
	
	//Check phonon supercell validity:
	if(wannier.phononSup.length_squared())
	{	//--- Check that supercell is diagonal in lattice directions:
		int superOffDiag = 0.;
		for(int j1=0; j1<3; j1++)
			for(int j2=0; j2<3; j2++)
				if(j1 != j2)
					superOffDiag += std::pow(supercell.super(j1,j2), 2);
		if(e.eInfo.qnums[0].k.length_squared() || superOffDiag)
			die("Phonon matrix elements require a Gamma-centered uniform kpoint mesh.\n");
		//--- Check that phonon supercell tiles the Wannier supercell:
		for(int j=0; j<3; j++)
		{	if(!wannier.phononSup[j] || supercell.super(j,j) % wannier.phononSup[j])
			{	die("Wannier supercell count %d is not a multiple of phonon supercell count %d for lattice direction %d.\n",
					supercell.super(j,j), wannier.phononSup[j], j);
			}
		}
		//--- Count phonon modes:
		nPhononModes = 0;
		invsqrtM.clear();
		string fname = wannier.getFilename(Wannier::FilenameInit, "phononBasis");
		if(mpiWorld->isHead())
		{	std::ifstream ifs(fname.c_str());
			if(ifs.is_open())
			{	while(!ifs.eof())
				{	string line; getline(ifs, line);
					trim(line);
					if(!line.length()) continue;
					istringstream iss(line);
					string spName; int atom; vector3<> disp; double M;
					iss >> spName >> atom >> disp[0] >> disp[1] >> disp[2] >> M;
					if(!iss.fail())
					{	nPhononModes++;
						invsqrtM.push_back(1./sqrt(M*amu));
					}
				}
			}
		}
		mpiWorld->bcast(nPhononModes);
		invsqrtM.resize(nPhononModes);
		mpiWorld->bcastData(invsqrtM);
		if(!nPhononModes) die("Error reading phonon modes from '%s'\n", fname.c_str());
		logPrintf("Found %d phonon modes in '%s'\n", nPhononModes, fname.c_str());
	}
}

void WannierMinimizer::initTransformDependent()
{
	//Initialize common basis:
	double kMaxSq = 0;
	for(const Kpoint& kpoint: kpoints)
		kMaxSq = std::max(kMaxSq, e.gInfo.GGT.metric_length_squared(kpoint.k));
	double GmaxEff = sqrt(2.*e.cntrl.Ecut) + sqrt(kMaxSq);
	double EcutEff = 0.5*GmaxEff*GmaxEff * (1.+symmThreshold); //add some margin for round-off error safety
	basis.setup(e.gInfo, e.iInfo, EcutEff, vector3<>());
	basisWrapper = std::make_shared<ColumnBundleTransform::BasisWrapper>(basis);
	
	//Initialize supercell basis (if necessary):
	if(needSuper)
	{	basisSuper.setup(gInfoSuper, e.iInfo, e.cntrl.Ecut, qnumSuper.k);
		basisSuperWrapper = std::make_shared<ColumnBundleTransform::BasisWrapper>(basisSuper);
	}
	
	//Initialize transforms:
	for(const Kpoint& kpoint: kpoints)
	{	const Basis& basisC = e.basis[kpoint.iReduced];
		const vector3<>& kC = e.eInfo.qnums[kpoint.iReduced].k;
		const matrix3<int>& super = e.coulombParams.supercell->super;
		transformMap[kpoint] = std::make_shared<ColumnBundleTransform>(kC, basisC, kpoint.k, *basisWrapper, nSpinor, sym[kpoint.iSym], kpoint.invert);
		if(needSuper)
			transformMapSuper[kpoint] = std::make_shared<ColumnBundleTransform>(kC, basisC, qnumSuper.k, *basisSuperWrapper, nSpinor, sym[kpoint.iSym], kpoint.invert, super);
	}
	
	//Read numerical trial orbitals if provided:
	if(wannier.numericalOrbitalsFilename.length())
	{	logPrintf("\n--- Initializing numerical trial orbitals ---\n");
		
		//Read the header:
		string fname = wannier.numericalOrbitalsFilename + ".header";
		std::ifstream ifs(fname.c_str());
		if(!ifs.is_open()) die("Could not open '%s' for reading.\n", fname.c_str());
		string comments;
		//--- Line 1:
		int nCols; size_t colLength;
		ifs >> nCols >> colLength;
		getline(ifs, comments);
		logPrintf("%d orbitals with %lu basis elements each.\n", nCols, colLength);
		//--- Line 2:
		matrix3<> GT;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ifs >> GT(i,j);
		getline(ifs, comments);
		//--- Basis elements:
		std::vector< vector3<int> > iGarr(colLength);
		vector3<int> hlfSin;
		for(vector3<int>& iG: iGarr)
		{	ifs >> iG[0] >> iG[1] >> iG[2];
			for(int k=0; k<3; k++) hlfSin[k] = std::max(hlfSin[k], abs(iG[k]));
		}
		ifs.close();
		
		//Create input grid:
		GridInfo gInfoIn;
		for(int k=0; k<3; k++)
		{	gInfoIn.S[k] = 2*(2*hlfSin[k] + 1);
			while(!fftSuitable(gInfoIn.S[k])) gInfoIn.S[k] += 2; 
		}
		gInfoIn.R = 2*M_PI * inv(~GT);
		gInfoIn.initialize(true);
		
		//Create input basis:
		Basis basisIn;
		std::vector<int> indexIn(colLength);
		for(size_t i=0; i<colLength; i++)
			indexIn[i] = gInfoIn.fullGindex(iGarr[i]);
		basisIn.setup(gInfoIn, e.iInfo, indexIn);
		
		//Read the wavefunctions:
		fname = wannier.numericalOrbitalsFilename;
		logPrintf("Reading numerical orbitals from '%s' ... ", fname.c_str()); logFlush();
		ColumnBundle Cin(nCols, colLength*nSpinor, &basisIn, &qnumSuper);
		Cin.read(fname.c_str());
		logPrintf("done.\n"); logFlush();
		
		//Convert wavefunctions to supercell basis:
		logPrintf("Converting numerical orbitals to supercell basis ... "); logFlush();
		//--- offset in input
		Cin = translate(Cin, -wannier.numericalOrbitalsOffset);
		//--- resample band-by-band"
		logSuspend();
		BlipResampler resample(gInfoIn, gInfoSuper);
		logResume();
		ColumnBundle C(nCols, basisSuper.nbasis*nSpinor, &basisSuper, &qnumSuper, isGpuEnabled());
		for(int b=0; b<nCols; b++)
			for(int s=0; s<nSpinor; s++)
				C.setColumn(b,s, J(resample(Cin.getColumn(b,s))));
		logPrintf("done.\n"); logFlush();
		
		//Split supercell wavefunction into kpoints:
		logPrintf("Dividing supercell numerical orbitals to k-points ... "); logFlush();
		for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,0) || isMine_q(ik,1))
		{	const KmeshEntry& ki = kMesh[ik];
			const ColumnBundleTransform& transform = *(transformMap.find(ki.point)->second);
			const ColumnBundleTransform& transformSuper = *(transformMapSuper.find(ki.point)->second);
			//Collect in k-dependent unit-cell basis:
			const Basis& basisC = e.basis[ki.point.iReduced];
			ColumnBundle temp(nCols, basisC.nbasis*nSpinor, &basisC, 0, isGpuEnabled());
			temp.zero();
			transformSuper.gatherAxpy(1., C,0,1, temp);
			//Transfer to common unit-cell basis:
			auto Ck = std::make_shared<ColumnBundle>(nCols, basis.nbasis*nSpinor, &basis, &ki.point, isGpuEnabled());
			Ck->zero();
			transform.scatterAxpy(1./ki.point.weight, temp, *Ck,0,1);
			numericalOrbitals[ki.point] = Ck;
		}
		logPrintf("done.\n"); logFlush();
	}
}

void WannierMinimizer::initRotations(int iSpin)
{
	//Load initial rotations if necessary:
	bool rotationsLoaded = false;
	if(wannier.loadRotations)
	{	//Read U
		string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU", &iSpin);
		bool Uexists = (fileSize(fname.c_str()) > 0);
		mpiWorld->bcast(Uexists); //Ensure MPI consistency of file check (avoid occassional NFS errors)
		if(Uexists)
		{	logPrintf("Reading initial rotations from '%s' ... ", fname.c_str()); logFlush();
			MPIUtil::File fp; mpiWorld->fopenRead(fp, fname.c_str());
			for(auto& ke: kMesh)
			{	ke.U.init(nBands, nCenters);
				mpiWorld->freadData(ke.U, fp);
			}
			mpiWorld->fclose(fp);
			logPrintf("done.\n"); logFlush();
			logPrintf("NOTE: ignoring trial orbitals since we are resuming a previous WannierMinimize.\n");
			rotationsLoaded = true;
		}
		else
		{	logPrintf("NOTE: no initial rotations in '%s'; using trial orbitals for first run.\n", fname.c_str());
			logFlush();
		}
	}
	
	//Load frozen rotations if necessary:
	std::vector<matrix> Ufrozen(kMesh.size());
	if(nFrozen)
	{	//Read U
		logPrintf("Reading frozen rotations from '%s' ... ", wannier.frozenUfilename.c_str()); logFlush();
		MPIUtil::File fp; mpiWorld->fopenRead(fp, wannier.frozenUfilename.c_str());
		for(size_t ik=0; ik<kMesh.size(); ik++)
		{	Ufrozen[ik].init(nBands, nFrozen);
			mpiWorld->freadData(Ufrozen[ik], fp);
		}
		mpiWorld->fread(rPinned.data(), sizeof(vector3<>), nFrozen, fp);
		mpiWorld->fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Compute / check the initial rotations:
	ostringstream ossErr;
	for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
	{	KmeshEntry& ke = kMesh[ik];
		string kString;
		{	ostringstream kOss;
			kOss << " at k = [ " << ke.point.k[0] << ' ' << ke.point.k[1] << ' ' << ke.point.k[2] << " ]\n";
			kString = kOss.str();
		}
		//Band ranges:
		int bStart=0, bStop=0, bFixedStart=0, bFixedStop=0;
		const std::vector<double>& eigs = e.eVars.Hsub_eigs[ke.point.iReduced + iSpin*qCount];
		if(wannier.outerWindow)
		{	bStart = 0;
			while(bStart<nBands && eigs[bStart]<wannier.eOuterMin)
				bStart++;
			bStop = bStart;
			while(bStop<nBands && eigs[bStop]<=wannier.eOuterMax)
				bStop++;
			if(bStop-bStart < nCenters)
			{	ossErr << "Number of bands within outer window = " << bStop-bStart
					<< " less than nCenters = " << nCenters << kString;
				break;
			}
			//Optionally range for inner window:
			if(wannier.innerWindow)
			{	bFixedStart = bStart;
				while(bFixedStart<bStop && eigs[bFixedStart]<wannier.eInnerMin)
					bFixedStart++;
				bFixedStop = bFixedStart;
				while(bFixedStop<bStop && eigs[bFixedStop]<=wannier.eInnerMax)
					bFixedStop++;
				if(bFixedStop-bFixedStart > nCenters)
				{	ossErr << "Number of bands within inner window = " << bFixedStop-bFixedStart
						<< " exceeds nCenters = " << nCenters << kString;
					break;
				}
			}
			else bFixedStart = bFixedStop = bStart; //fixed interval is empty
		}
		else //fixed bands
		{	bFixedStart = bStart = wannier.bStart;
			bFixedStop  = bStop  = wannier.bStart + nCenters;
		}
		
		//Check compatibility of frozen window with band range:
		if(nFrozen)
		{	const double tol = 1e-6 * nFrozen;
			if( (bFixedStart>0 && nrm2(Ufrozen[ik](0,bFixedStart, 0,nFrozen)) > tol)
			 || (bFixedStop<nBands && nrm2(Ufrozen[ik](bFixedStop,nBands, 0,nFrozen)) > tol) )
			{	ossErr << "Frozen rotations are incompatible with current outer window" << kString;
				break;
			}
		}
		
		//Initial rotation of bands to get to Wannier subspace:
		ke.nFixed = bFixedStop - bFixedStart;
		int nFree = nCenters - ke.nFixed;
		ke.nIn = (nFree>0) ? (bStop-bStart) : nCenters; //number of bands contributing to Wannier subspace
		const double tol = 1e-6 * nCenters; //tolerance for checking unitarity etc.
		//--- Create fixed part of U1 (if any):
		ke.U1 = zeroes(nBands, ke.nIn);
		if(ke.nFixed > 0)
		{	matrix U1fixed = eye(ke.nFixed);
			if(nFrozen)
			{	if(ke.nFixed < nFrozen)
				{	ossErr << "Number of bands within " << (wannier.innerWindow ? "inner window" : "fixed band range")
						<< " = " << bStop-bStart << " less than nFrozen = " << nFrozen << kString;
					break;
				}
				//Rotate the fixed subspace so that the frozen centers come first:
				matrix UfrozenNZ = Ufrozen[ik](bFixedStart,bFixedStop, 0,nFrozen); //drop the zero rows (should be zero, else SVD check below will fail)
				U1fixed.set(0,ke.nFixed, 0,nFrozen, UfrozenNZ);
				//fill in the null space of UfrozenNZ using an SVD
				matrix U, Vdag; diagMatrix S;
				UfrozenNZ.svd(U, S, Vdag);
				if(nrm2(S(0,nFrozen)-eye(nFrozen))>tol)
				{	ossErr << "Frozen rotations are incompatible with current outer window" << kString;
					break;
				}
				if(nFrozen<ke.nFixed)
					U1fixed.set(0,ke.nFixed, nFrozen,ke.nFixed, U(0,ke.nFixed, nFrozen,ke.nFixed));
			}
			ke.U1.set(bFixedStart,bFixedStop, 0,ke.nFixed, U1fixed);
		}
		if(rotationsLoaded)
		{	//Factorize U (nBands x nCenters) into U1 (nBands x nIn) and U2 (nIn x nIn):
			//--- check unitarity:
			if(nrm2(dagger(ke.U) * ke.U - eye(nCenters)) > tol)
			{	ossErr << "Initial rotations U are not unitary" << kString;
				break;
			}
			//--- check compatibility with frozen rotations:
			if(nFrozen && nrm2(ke.U(0,nBands, 0,nFrozen) - Ufrozen[ik])>tol)
			{	ossErr << "Initial rotations U are incompatible with frozen rotations" << kString;
				break;
			}
			//--- determine the free columns of U1:
			if(nFree > 0)
			{	matrix U1free = ke.U(bStart,bStop, 0,nCenters);
				//project out the fixed columns (if any):
				if(ke.nFixed > 0)
				{	matrix U1fixed = ke.U1(bStart,bStop, 0,ke.nFixed);
					U1free -= U1fixed*(dagger(U1fixed) * U1free);
					//SVD to find the remaining non-null columns:
					matrix U, Vdag; diagMatrix S;
					U1free.svd(U, S, Vdag);
					if(nrm2(S(0,nFree)-eye(nFree))>tol || nrm2(S(nFree,nCenters))>tol)
					{	ossErr << "Initial rotations are incompatible with current inner window" << kString;
						break;
					}
					U1free = U(0,ke.nIn, 0,nFree); //discard null subspace
				}
				ke.U1.set(bStart,bStop, ke.nFixed,nCenters, U1free);
			}
			//--- check U1:
			if(nrm2((dagger(ke.U1) * ke.U1)(0,nCenters, 0,nCenters) - eye(nCenters)) > tol)
			{	ossErr << "Initial rotations U1 are not unitary" << kString;
				break;
			}
			if( (bStart>0 && nrm2(ke.U1(0,bStart, 0,nCenters))>tol) 
				|| (bStop<nBands && nrm2(ke.U1(bStop,nBands, 0,nCenters))>tol) )
			{	ossErr << "Initial rotations are incompatible with current outer window / band selection" << kString;
				break;
			}
			//--- calculate and check U2:
			ke.U2 = eye(ke.nIn);
			ke.U2.set(nFrozen,nCenters, nFrozen,nCenters,
				dagger(ke.U1(0,nBands, nFrozen,nCenters)) * ke.U(0,nBands, nFrozen,nCenters));
			if(nrm2(dagger(ke.U2) * ke.U2 - eye(ke.nIn)) > tol)
			{	ossErr << "Initial rotations U2 are not unitary" << kString;
				break;
			}
			//--- Compute extra linearly indep columns of U1 (if necessary):
			if(ke.nIn > nCenters)
			{	matrix U, Vdag; diagMatrix S;
				matrix(ke.U1(bStart,bStop, 0,ke.nIn)).svd(U, S, Vdag);
				ke.U1.set(bStart,bStop, 0,ke.nIn, U*Vdag);
			}
		}
		else
		{	//Determine from trial orbitals:
			matrix CdagG = getWfns(ke.point, iSpin) ^ O(trialWfns(ke.point));
			int nNew = nCenters - nFrozen; //number of new centers
			//--- Pick up best linear combination of remaining bands (if any)
			if(nFree > 0)
			{	//Create overlap matrix with contribution from fixed bands projected out:
				matrix CdagGFree;
				if(ke.nFixed > 0)
				{	//SVD the fixed band contribution to the trial space:
					matrix U, Vdag; diagMatrix S;
					matrix(CdagG(bFixedStart,bFixedStop, 0,nNew)).svd(U, S, Vdag);
					//Project out the fixed bands (use only the zero singular values)
					CdagGFree = CdagG * dagger(Vdag(nNew-nFree,nNew, 0,nNew));
				}
				else CdagGFree = CdagG;
				//Truncate to non-zero rows:
				int nLo = bFixedStart-bStart;
				int nHi = bStop-bFixedStop;
				int nOuter = nLo+nHi;
				matrix CdagGFreeNZ = zeroes(nOuter, nOuter);
				if(nLo>0) CdagGFreeNZ.set(0,nLo, 0,nFree, CdagGFree(bStart,bFixedStart, 0,nFree));
				if(nHi>0) CdagGFreeNZ.set(nLo,nOuter, 0,nFree, CdagGFree(bFixedStop,bStop, 0,nFree));
				//SVD to get best linear combinations first:
				matrix U, Vdag; diagMatrix S;
				CdagGFreeNZ.svd(U, S, Vdag);
				//Convert left space from non-zero to all bands:
				matrix Upad = zeroes(nBands, nOuter);
				if(nLo>0) Upad.set(bStart,bFixedStart, 0,nOuter, U(0,nLo, 0,nOuter));
				if(nHi>0) Upad.set(bFixedStop,bStop, 0,nOuter, U(nLo,nOuter, 0,nOuter));
				//Include this combination in U1:
				ke.U1.set(0,nBands, ke.nFixed,ke.nIn, Upad * Vdag);
			}
			
			//Optimal initial rotation within Wannier subspace:
			matrix WdagG = dagger(ke.U1(0,nBands, nFrozen,nCenters)) * CdagG;
			ke.U2 = eye(ke.nIn);
			bool isSingular = false;
			ke.U2.set(nFrozen,nCenters, nFrozen,nCenters, WdagG * invsqrt(dagger(WdagG) * WdagG, 0, 0, &isSingular));
			if(isSingular)
			{	ossErr << "Trial orbitals are linearly dependent / do not cover desired subspace" << kString;
				break;
			}
		}
		//Make initial rotations exactly unitary:
		bool isSingular = false;
		ke.U1 = fixUnitary(ke.U1, &isSingular);
		ke.U2 = fixUnitary(ke.U2, &isSingular);
		if(isSingular)
		{	ossErr << "Initial rotations are singular" << kString;
			break;
		}
		ke.U = ke.U1 * ke.U2;
	}
	mpiWorld->checkErrors(ossErr);
	suspendOperatorThreading();
	
	//Broadcast initial rotations:
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	KmeshEntry& ke = kMesh[ik];
		mpiWorld->bcast(ke.nIn, whose_q(ik,iSpin));
		mpiWorld->bcast(ke.nFixed, whose_q(ik,iSpin));
		if(!isMine_q(ik,iSpin))
		{	ke.U1 = zeroes(nBands, ke.nIn);
			ke.U2 = zeroes(ke.nIn, ke.nIn); //truncated to nCenters x nCenters after minimization
			ke.U = zeroes(nBands, ke.nIn); //truncated to nBands x nCenters after minimization
		}
		mpiWorld->bcastData(ke.U1, whose_q(ik,iSpin));
		mpiWorld->bcastData(ke.U2, whose_q(ik,iSpin));
		mpiWorld->bcastData(ke.U, whose_q(ik,iSpin));
		if(!isMine(ik)) //No longer need sub-matrices on this process
		{	ke.U1 = matrix();
			ke.U2 = matrix();
			if(!ke.mpi) ke.U = matrix(); //No longer need U on this process either
		}
	}
}
