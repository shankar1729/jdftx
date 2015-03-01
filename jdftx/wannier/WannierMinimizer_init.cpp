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
#include <electronic/Blip.h>
#include <electronic/operators.h>

WannierMinimizer::WannierMinimizer(const Everything& e, const Wannier& wannier, bool needSuperOverride)
: e(e), wannier(wannier), sym(e.symm.getMatrices()),
	nCenters(wannier.nCenters), nFrozen(wannier.nFrozen), nBands(e.eInfo.nBands),
	nSpins(e.eInfo.spinType==SpinZ ? 2 : 1), qCount(e.eInfo.qnums.size()/nSpins),
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

	//Initialize cell map (for matrix element output):
	iCellMap = getCellMap(e.gInfo.R, gInfoSuper.R, wannier.getFilename(Wannier::FilenameDump, "mlwfCellMap"));
	
	//Determine and output the band ranges:
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{	std::vector<double> eMin(nBands, DBL_MAX), eMax(nBands, -DBL_MAX);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) if(e.eInfo.qnums[q].spin==iSpin)
			for(int b=0; b<nBands; b++)
			{	eMin[b] = std::min(eMin[b], e.eVars.Hsub_eigs[q][b]);
				eMax[b] = std::max(eMax[b], e.eVars.Hsub_eigs[q][b]);
			}
		mpiUtil->allReduce(eMin.data(), eMin.size(), MPIUtil::ReduceMin);
		mpiUtil->allReduce(eMax.data(), eMax.size(), MPIUtil::ReduceMax);
		if(mpiUtil->isHead())
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
		ki.point.k = (~sym[src.iSym]) * (src.invert * qnums[src.iReduced].k) + src.offset;
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
	kDivision.init(kMesh.size(), mpiUtil);
	kDivision.myRange(ikStart, ikStop);
	
	//Initialized pinned arrays:
	for(int n=nFrozen; n<nCenters; n++)
	{	pinned[n] = wannier.trialOrbitals[n-nFrozen].pinned;
		if(pinned[n])
			rPinned[n] = e.gInfo.R * wannier.trialOrbitals[n-nFrozen].rCenter;
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
		//--- Generate phonon cell map and output cellMapSq:
		phononCellMap = getCellMap(e.gInfo.R, e.gInfo.R * Diag(wannier.phononSup));
		if(mpiUtil->isHead())
		{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellMapSqPh");
			logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
			FILE* fp = fopen(fname.c_str(), "w");
			fprintf(fp, "#i0 i1 i2  i0' i1' i2'   (integer lattice combinations for pairs of sites)\n");
			for(const auto& entry1: phononCellMap)
			for(const auto& entry2: phononCellMap)
			{	const vector3<int>& iR1 = entry1.first;
				const vector3<int>& iR2 = entry2.first;
				fprintf(fp, "%+2d %+2d %+2d  %+2d %+2d %+2d\n", iR1[0], iR1[1], iR1[2], iR2[0], iR2[1], iR2[2]);
			}
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
		//--- Count phonon modes:
		nPhononModes = 0;
		string fname = wannier.getFilename(Wannier::FilenameInit, "phononBasis");
		if(mpiUtil->isHead())
		{	std::ifstream ifs(fname.c_str());
			if(ifs.is_open())
			{	while(!ifs.eof())
				{	string line; getline(ifs, line);
					trim(line);
					if(!line.length()) continue;
					istringstream iss(line);
					string spName; int atom; vector3<> disp;
					iss >> spName >> atom >> disp[0] >> disp[1] >> disp[2];
					if(!iss.fail()) nPhononModes++;
				}
			}
		}
		mpiUtil->bcast(nPhononModes);
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
		{	ColumnBundle temp(nCols, basis.nbasis, &basis, 0, isGpuEnabled());
			for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,0) || isMine_q(ik,1))
			{	const KmeshEntry& ki = kMesh[ik];
				const ColumnBundleTransform& transform = *(transformMap.find(ki.point)->second);
				const ColumnBundleTransform& transformSuper = *(transformMapSuper.find(ki.point)->second);
				auto Ck = std::make_shared<ColumnBundle>(nCols, basis.nbasis*nSpinor, &basis, &ki.point, isGpuEnabled());
				temp.zero();
				transformSuper.gatherAxpy(1., C,0,1, temp);
				Ck->zero();
				transform.scatterAxpy(1./ki.point.weight, temp, *Ck,0,1);
				numericalOrbitals[ki.point] = Ck;
			}
		}
		logPrintf("done.\n"); logFlush();
	}
}
