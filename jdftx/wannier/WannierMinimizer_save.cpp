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

#include <wannier/WannierMinimizerFD.h>
#include <core/ScalarFieldIO.h>
#include <core/SphericalHarmonics.h>

void WannierMinimizer::saveMLWF()
{	for(int iSpin: wannier.iSpinArr)
	{	logPrintf("\n");
		saveMLWF(iSpin);
	}
}

void WannierMinimizer::saveMLWF(int iSpin)
{	
	//Load / compute and broadcast initial rotations:
	initRotations(iSpin);
	
	//Sub-class specific initialization:
	initialize(iSpin);
	
	//Minimize:
	double Omega = minimize(wannier.minParams);
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	KmeshEntry& ki = kMesh[ik];
		//Make U available on all processes:
		if(isMine(ik))
			ki.U = ki.U(0,nBands, 0,nCenters); //no longer need unitary completion after minimize
		else
			ki.U = zeroes(nBands, nCenters);
		mpiWorld->bcastData(ki.U, whose(ik));
		if(ki.U2) ki.U2 = ki.U2(0,nCenters, 0,nCenters); //no longer need unitary completion after minimize
	}
	double OmegaI = getOmegaI();
	logPrintf("\nOptimum spread:\n\tOmega:  %.15le\n\tOmegaI: %.15le\n", Omega, OmegaI);
	
	//List the centers:
	logPrintf("\nCenters in %s coords:\n", e.iInfo.coordsType==CoordsCartesian ? "cartesian" : "lattice");
	for(int n=0; n<nCenters; n++)
	{	vector3<> rCoords = e.iInfo.coordsType==CoordsCartesian
			? rExpect[n] : e.gInfo.invR * rExpect[n]; //r in coordinate system of choice
		logPrintf("\t[ %lg %lg %lg ] spread: %lg bohr^2\n", rCoords[0], rCoords[1], rCoords[2], rSqExpect[n] - rExpect[n].length_squared());
	}
	logFlush();
	
	//Save the matrices:
	if(mpiWorld->isHead() && wannier.minParams.nIterations) //re-save only if any minimization has occured
	{	//Write U:
		string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str());
		FILE* fp = fopen(fname.c_str(), "w");
		for(const auto& ke: kMesh) ke.U.write(fp);
		fwrite(rPinned.data(), sizeof(vector3<>), nCenters, fp); //read back only if used as frozen centers
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output range of centers that span each band
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfBandContrib", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		//Find range of centers that span each band:
		std::vector<int> nMin(nBands, nCenters), nMax(nBands, 0);
		for(const KmeshEntry& ke: kMesh)
			for(int n=0; n<nCenters; n++)
				for(int b=0; b<nBands; b++)
					if(ke.U(b,n).norm() > 1e-12)
					{	nMin[b] = std::min(nMin[b], n);
						nMax[b] = std::max(nMax[b], n);
					}
		//Find energy range that is exact for each range of centers:
		std::map< std::pair<int,int>, std::pair<double,double> > nRangeToErange;
		nRangeToErange[std::make_pair(-1,-1)] = std::make_pair(NAN,NAN);
		for(int b=0; b<nBands; b++)
			if(nMin[b] <= nMax[b])
			{	std::pair<int,int> nRange(nMin[b], nMax[b]);
				if(nRangeToErange.find(nRange) == nRangeToErange.end())
				{	double eMin =-INFINITY, eMax = +INFINITY;
					for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
					{	KmeshEntry& ke = kMesh[ik];
						const std::vector<double>& eigs = e.eVars.Hsub_eigs[ke.point.iReduced + iSpin*qCount];
						matrix Utrunc = ke.U(0,nBands, nMin[b],nMax[b]+1);
						diagMatrix overlap = diag(Utrunc * dagger(Utrunc));
						double eMin_k = +INFINITY, eMax_k = -INFINITY;
						for(int iBand=0; iBand<nBands; iBand++)
							if(fabs(overlap[iBand] - 1.) < 1e-6)
							{	eMin_k = std::min(eMin_k, iBand ? eigs[iBand-1] : -INFINITY);
								eMax_k = std::max(eMax_k, (iBand<nBands) ? eigs[iBand+1] : +INFINITY);
							}
						eMin = std::max(eMin, eMin_k);
						eMax = std::min(eMax, eMax_k);
					}
					mpiWorld->allReduce(eMin, MPIUtil::ReduceMax);
					mpiWorld->allReduce(eMax, MPIUtil::ReduceMin);
					if(eMin >= eMax) { eMin = eMax = NAN; }
					nRangeToErange[nRange] = std::make_pair(eMin, eMax);
				}
			}
			else { nMin[b] = nMax[b] = -1; } //unused band
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(int b=0; b<nBands; b++)
			{	const std::pair<double,double> eRange = nRangeToErange[std::make_pair(nMin[b],nMax[b])];
				fprintf(fp, "%d %d   %+10.5lf %+10.5lf\n", nMin[b], nMax[b], eRange.first, eRange.second);
			}
			fclose(fp);
		}
		logPrintf("done.\n"); logFlush();
	}
	
	realPartOnly = !e.eInfo.isNoncollinear(); //cannot save only real part in noncollinear calculations
	
	//Save wannier wave functions if needed:
	if(wannier.saveWfns || wannier.saveWfnsRealSpace)
		saveMLWF_C(iSpin);
	
	//Initialize cell map (for matrix element output):
	//--- get wannier centers in lattice coordinates:
	xExpect.clear();
	for(vector3<> r: rExpect)
		xExpect.push_back(e.gInfo.invR * r);
	//--- get cell map with weights based on these center positions:
	iCellMap = getCellMap(e.gInfo.R, gInfoSuper.R, e.coulombParams.isTruncated(),
		xExpect, xExpect, wannier.rSmooth, wannier.getFilename(Wannier::FilenameDump, "mlwfCellMap", &iSpin));
	//--- output cell map weights:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellWeights", &iSpin);
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("could not open file for writing.\n");
		for(auto iter: iCellMap)
			iter.second.write_real(fp);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Create unique cells
	vector3<int> kfold = e.eInfo.kFoldingCount();
	int kfoldProd = kfold[0]*kfold[1]*kfold[2];
	std::vector<vector3<int>> uniqueCells; uniqueCells.reserve(kfoldProd);
	{	vector3<int> iR;
		for(iR[0]=0; iR[0]<kfold[0]; iR[0]++)
		for(iR[1]=0; iR[1]<kfold[1]; iR[1]++)
		for(iR[2]=0; iR[2]<kfold[2]; iR[2]++)
			uniqueCells.push_back(iR);
	}
	
	//Compute Fourier transform phase (common to all electronic k-mesh outputs below):
	int nqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		nqMine++;
	nqMine = std::max(nqMine,1); //avoid zero size matrices below
	matrix phase = zeroes(nqMine, uniqueCells.size());
	{	int iqMine = 0;
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	int iCell = 0;
			for(const vector3<int>& iR: uniqueCells)
				phase.set(iqMine, iCell++, kMesh[i].point.weight * cis(-2*M_PI*dot(kMesh[i].point.k, iR)));
			iqMine++;
		}
	}
	resumeOperatorThreading();
	
	//Electronic k-mesh outputs:
	saveMLWF_H(iSpin, phase); //Hamiltonian
	if(wannier.saveMomenta) saveMLWF_P(iSpin, phase); //Momenta
	if(wannier.saveSpin) saveMLWF_S(iSpin, phase); //Spins
	if(wannier.saveR or wannier.saveRP) saveMLWF_RP(iSpin, phase); //R and/or R*P matrix elements
	if(wannier.zVfilename.length()) saveMLWF_Z(iSpin, phase); //z position
	if(wannier.zH) saveMLWF_W(iSpin, phase); //Slab weights
	saveMLWF_ImSigma_ee(iSpin, phase);
	
	//Phonon q-mesh outputs:
	bool savePhonon = wannier.phononSup.length_squared();
	if(savePhonon) saveMLWF_D(iSpin, phase); //Gradient (for phonon sum rule)
	if(savePhonon) saveMLWF_phonon(iSpin);
	
	//Defect outputs:
	for(const DefectSupercell& ds: wannier.defects)
		saveMLWF_defect(iSpin, (DefectSupercell&)ds);
	suspendOperatorThreading();
}


//-------- Function implementations for outputting each supported quantity above ---------


void WannierMinimizer::saveMLWF_C(int iSpin)
{	resumeOperatorThreading();
	
	//Compute supercell wavefunctions:
	logPrintf("Computing supercell wavefunctions ... "); logFlush();
	ColumnBundle Csuper(nCenters, basisSuper.nbasis*nSpinor, &basisSuper, &qnumSuper, isGpuEnabled());
	Csuper.zero();
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	const KmeshEntry& ki = kMesh[i];
		axpyWfns(ki.point.weight, ki.U, ki.point, iSpin, Csuper);
	}
	mpiWorld->allReduceData(Csuper, MPIUtil::ReduceSum);
	Csuper = translate(Csuper, vector3<>(.5,.5,.5)); //center in supercell
	logPrintf("done.\n"); logFlush();
	
	//Save supercell wavefunctions in reciprocal space:
	if(mpiWorld->isHead() && wannier.saveWfns)
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfC", &iSpin);
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
		Csuper.write(fname.c_str());
		logPrintf("done.\n"); logFlush();
		//Header:
		fname = fname + ".header";
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "%d %lu #nColumns, columnLength\n", nCenters, basisSuper.nbasis);
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				fprintf(fp, "%.15g ", gInfoSuper.GT(i,j));
		fprintf(fp, "#GT row-major (G col-major), iGarr follows:\n");
		for(const vector3<int>& iG: basisSuper.iGarr)
			fprintf(fp, "%d %d %d\n", iG[0], iG[1], iG[2]);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Save supercell wavefunctions in real space:
	if(mpiWorld->isHead() && wannier.saveWfnsRealSpace) for(int n=0; n<nCenters; n++) for(int s=0; s<nSpinor; s++)
	{	//Generate filename
		ostringstream varName;
		varName << (nSpinor*n+s) << ".mlwf";
		string fname = wannier.getFilename(Wannier::FilenameDump, varName.str(), &iSpin);
		logPrintf("Dumping '%s':\n", fname.c_str());
		//Convert to real space and optionally remove phase:
		complexScalarField psi = I(Csuper.getColumn(n,s));
		if(qnumSuper.k.length_squared() > symmThresholdSq)
			multiplyBlochPhase(psi, qnumSuper.k);
		if(realPartOnly)
		{	complex* psiData = psi->data();
			double meanPhase, sigmaPhase, rmsImagErr;
			removePhase(gInfoSuper.nr, psiData, meanPhase, sigmaPhase, rmsImagErr);
			logPrintf("\tPhase = %lf +/- %lf\n", meanPhase, sigmaPhase); logFlush();
			logPrintf("\tRMS imaginary part = %le (after phase removal)\n", rmsImagErr);
			logFlush();
			//Write real part of supercell wavefunction to file:
			FILE* fp = fopen(fname.c_str(), "wb");
			if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
			for(int i=0; i<gInfoSuper.nr; i++)
				fwriteLE(psiData+i, sizeof(double), 1, fp);
			fclose(fp);
		}
		else saveRawBinary(psi, fname.c_str());
	}
	
	suspendOperatorThreading();
}


//Save Hamiltonian in Wannier basis:
void WannierMinimizer::saveMLWF_H(int iSpin, const matrix& phase)
{	matrix HwannierTilde = zeroes(nCenters*nCenters, phase.nRows());
	int iqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	//Apply MLWF-optimized unitary rotations to Hamiltonian at current k:
		matrix Hsub = dagger(kMesh[i].U) * e.eVars.Hsub_eigs[kMesh[i].point.iReduced + iSpin*qCount] * kMesh[i].U;
		callPref(eblas_copy)(HwannierTilde.dataPref()+HwannierTilde.index(0,iqMine), Hsub.dataPref(), Hsub.nData());
		iqMine++;
	}
	//Fourier transform to Wannier space and save
	dumpWannierized(HwannierTilde, phase, "mlwfH", realPartOnly, iSpin);
}


//Save momenta in Wannier basis:
void WannierMinimizer::saveMLWF_P(int iSpin, const matrix& phase)
{	assert(wannier.saveMomenta);
	//Compute momentum matrix elements of Bloch states:
	std::vector<vector3<matrix>> pBloch(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		if(e.eInfo.qnums[q].index()==iSpin)
			for(int iDir=0; iDir<3; iDir++)
				pBloch[q][iDir] = e.iInfo.rHcommutator(e.eVars.C[q], iDir, e.eVars.Hsub_eigs[q]); //note factor of -iota dropped to make it real (and anti-symmetric)
	//Convert to Wannier basis:
	matrix pWannierTilde = zeroes(nCenters*nCenters*3, phase.nRows());
	int iqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	matrix pSub(nCenters*nCenters, 3);
		for(int iDir=0; iDir<3; iDir++)
		{	matrix pSubDir = pBloch[kMesh[i].point.iReduced + iSpin*qCount][iDir];
			if(kMesh[i].point.invert<0) //apply complex conjugate:
				callPref(eblas_dscal)(pSubDir.nData(), -1., ((double*)pSubDir.dataPref())+1, 2);
			pSubDir = dagger(kMesh[i].U) * pSubDir * kMesh[i].U; //apply MLWF-optimized rotations
			callPref(eblas_copy)(pSub.dataPref()+pSub.index(0,iDir), pSubDir.dataPref(), pSubDir.nData());
		}
		//Store with spatial transformation:
		matrix rot(e.gInfo.R * sym[kMesh[i].point.iSym].rot * e.gInfo.invR); //cartesian symmetry matrix
		pSub = pSub * rot;
		callPref(eblas_copy)(pWannierTilde.dataPref()+pWannierTilde.index(0,iqMine), pSub.dataPref(), pSub.nData());
		iqMine++;
	}
	//Fourier transform to Wannier space and save
	dumpWannierized(pWannierTilde, phase, "mlwfP", realPartOnly, iSpin);
}

//Save gradient matrix elements in Wannier basis:
void WannierMinimizer::saveMLWF_D(int iSpin, const matrix& phase)
{	bool savePhonon = wannier.phononSup.length_squared();
	assert(savePhonon);
	DblochMesh.resize(kMesh.size());
	//Compute momentum matrix elements of Bloch states:
	std::vector<vector3<matrix>> Dbloch(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		if(e.eInfo.qnums[q].index()==iSpin)
			for(int iDir=0; iDir<3; iDir++)
				Dbloch[q][iDir] = e.gInfo.detR * (e.eVars.C[q] ^ D(e.eVars.C[q], iDir));
	//Convert to Wannier basis:
	matrix DwannierTilde = zeroes(nCenters*nCenters*3, phase.nRows());
	int iqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	matrix Dsub(nCenters*nCenters, 3);
		DblochMesh[i].init(nBands*nBands, 3);
		for(int iDir=0; iDir<3; iDir++)
		{	matrix DsubDir = Dbloch[kMesh[i].point.iReduced + iSpin*qCount][iDir];
			if(kMesh[i].point.invert<0) //apply complex conjugate:
				callPref(eblas_dscal)(DsubDir.nData(), -1., ((double*)DsubDir.dataPref())+1, 2);
			callPref(eblas_copy)(DblochMesh[i].dataPref()+DblochMesh[i].index(0,iDir), DsubDir.dataPref(), DsubDir.nData());  //copy before rotations
			DsubDir = dagger(kMesh[i].U) * DsubDir * kMesh[i].U; //apply MLWF-optimized rotations
			callPref(eblas_copy)(Dsub.dataPref()+Dsub.index(0,iDir), DsubDir.dataPref(), DsubDir.nData());
		}
		//Store with spatial transformation:
		matrix rot(e.gInfo.R * sym[kMesh[i].point.iSym].rot * e.gInfo.invR); //cartesian symmetry matrix
		DblochMesh[i] = DblochMesh[i] * rot;
		Dsub = Dsub * rot;
		callPref(eblas_copy)(DwannierTilde.dataPref()+DwannierTilde.index(0,iqMine), Dsub.dataPref(), Dsub.nData());
		iqMine++;
	}
	//Fourier transform to Wannier space and save
	dumpWannierized(DwannierTilde, phase, "mlwfD", realPartOnly, iSpin);
}


//Save spin in Wannier basis:
void WannierMinimizer::saveMLWF_S(int iSpin, const matrix& phase)
{	assert(wannier.saveSpin);
	assert(e.eInfo.isNoncollinear());
	//Compute spin matrix elements of Bloch states:
	std::vector<vector3<matrix>> Sbloch(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		if(e.eInfo.qnums[q].index()==iSpin)
			Sbloch[q] = spinOverlap(e.eVars.C[q]);
	//Convert to Wannier basis:
	matrix SwannierTilde = zeroes(nCenters*nCenters*3, phase.nRows());
	int iqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	matrix Ssub(nCenters*nCenters, 3);
		for(int iDir=0; iDir<3; iDir++)
		{	matrix SsubDir = Sbloch[kMesh[i].point.iReduced + iSpin*qCount][iDir];
			if(kMesh[i].point.invert<0) //apply negative complex conjugate (because spin is a pseudo-vector):
				callPref(eblas_dscal)(SsubDir.nData(), -1., ((double*)SsubDir.dataPref())+0, 2);
			SsubDir = dagger(kMesh[i].U) * SsubDir * kMesh[i].U; //apply MLWF-optimized rotations
			callPref(eblas_copy)(Ssub.dataPref()+Ssub.index(0,iDir), SsubDir.dataPref(), SsubDir.nData());
		}
		//Store with spatial transformation:
		matrix3<> rot = e.gInfo.R * sym[kMesh[i].point.iSym].rot * e.gInfo.invR; //cartesian symmetry matrix
		Ssub = Ssub * matrix(rot * (1./det(rot))); //extra rotation sign because spin is a pseudo-vector
		callPref(eblas_copy)(SwannierTilde.dataPref()+SwannierTilde.index(0,iqMine), Ssub.dataPref(), Ssub.nData());
		iqMine++;
	}
	//Fourier transform to Wannier space and save
	dumpWannierized(SwannierTilde, phase, "mlwfS", realPartOnly, iSpin);
}

//Save R and/or R*P matrix elements:
void WannierMinimizer::saveMLWF_RP(int iSpin, const matrix& phase)
{	assert(wannier.saveR or wannier.saveRP);

	//Compute r and/or r*p matrix elements from dC/dk on Wannier mesh:
	logPrintf("Computing <du/dk|...|u> on Wannier mesh ... "); logFlush();
	std::vector<vector3<matrix>> rMesh(kMesh.size());
	std::vector<matrix3<matrix>> rrHmesh(kMesh.size());
	assert(wannier.localizationMeasure == Wannier::LM_FiniteDifference);
	const WannierMinimizerFD& wmin_fd = *((const WannierMinimizerFD*)this);
	Cother.assign(e.eInfo.nStates, ColumnBundle());
	int qInterval = std::max(1, int(round(qCount/20.))); //interval for reporting progress
	for(int jReduced=0; jReduced<qCount; jReduced++)
	{	int qj = jReduced + qCount*iSpin; //current reduced state whose contributions to dC/dk at other states is being calculated
		//Get this wavefunction to all processes:
		int jProcess = e.eInfo.whose(qj);
		if(jProcess == mpiWorld->iProcess()) //send
			mpiWorld->bcastData((ColumnBundle&)e.eVars.C[qj], jProcess);
		else //recv
		{	Cother[qj].init(nBands, e.basis[qj].nbasis*nSpinor, &e.basis[qj], &e.eInfo.qnums[qj]);
			mpiWorld->bcastData(Cother[qj], jProcess);
		}
		//Find all edges connecting this one to states I have data for:
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	for(const WannierMinimizerFD::Edge& edge: wmin_fd.edges_bi[i])
				if(int(edge.point.iReduced) == jReduced)
				{	const KmeshEntry& ki = kMesh[i];
					const KmeshEntry& kj = kMesh[edge.ik];
					ColumnBundle Ci = getWfns(ki.point, iSpin) * ki.U;
					ColumnBundle Cj = getWfns(edge.point, iSpin) * kj.U;
					Cj.qnum = Ci.qnum; //effectively transport coefficients kj -> ki
					//Compute Cj ^ [r,H] Ci:
					int qi = ki.point.iReduced + qCount*iSpin;
					matrix Cj_O_Ci = Cj ^ O(Ci);
					vector3<matrix> Cj_rHcomm_Ci;
					if(wannier.saveRP)
					{	matrix Cj_H_Ci = Cj_O_Ci * (dagger(ki.U) * e.eVars.Hsub_eigs[qi] * ki.U);
						Cj_rHcomm_Ci = e.iInfo.rHcommutator(Cj, Ci, Cj_H_Ci);
					}
					//Collect contributions to (i dC/dk) ^ C  and/or (i dC/dk) ^ [r,H] C
					for(int mu=0; mu<3; mu++) //d/dk i.e. r direction
						if(edge.b[mu])
						{	complex w_mu(0, -0.5 * edge.wb * edge.b[mu]);
							if(wannier.saveR)
								rMesh[i][mu] += w_mu * Cj_O_Ci;
							if(wannier.saveRP)
								for(int nu=0; nu<3; nu++) //[r,H] i.e. p direction
									rrHmesh[i](mu, nu) += w_mu * Cj_rHcomm_Ci[nu];
						}
				}
		}
		Cother[qj] = 0; //free
		//Print progress:
		if((jReduced+1) % qInterval == 0)
		{	logPrintf("%d%% ", int(round((jReduced+1)*100./qCount)));
			logFlush();
		};
	}
	Cother.clear();
	logPrintf("done.\n"); logFlush();
	
	if(wannier.saveR)
	{	//Collect local contributions in contiguous form:
		int nCentersSq = nCenters * nCenters;
		matrix RwannierTilde = zeroes(nCentersSq*3, phase.nRows());
		complex* Rdata = RwannierTilde.dataPref();
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	for(int iDir=0; iDir<3; iDir++)
			{	callPref(eblas_copy)(Rdata, rMesh[i][iDir].dataPref(), nCentersSq);
				Rdata += nCentersSq;
			}
		}
		
		//Fourier transform to Wannier space and save:
		dumpWannierized(RwannierTilde, phase, "mlwfR", realPartOnly, iSpin);
	}
	
	if(wannier.saveRP)
	{	//Collect local contributions in contiguous form:
		int nCentersSq = nCenters * nCenters;
		matrix RPwannierTilde = zeroes(nCentersSq*3*3, phase.nRows());
		complex* RPdata = RPwannierTilde.dataPref();
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	for(int iDir=0; iDir<3; iDir++)
				for(int jDir=0; jDir<3; jDir++)
				{	callPref(eblas_copy)(RPdata, rrHmesh[i](iDir, jDir).dataPref(), nCentersSq);
					RPdata += nCentersSq;
				}
		}
		
		//Fourier transform to Wannier space and save:
		dumpWannierized(RPwannierTilde, phase, "mlwfRP", realPartOnly, iSpin);
	}
}


//Thread function for setting fourier transform of z-slice of half-width zH centered at z0 with smoothness zSigma:
inline void slabWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
	complex* w, double z0, double zH, double zSigma)
{	double prefac = 2*zH; //= average value of weight in unit cell
	THREAD_halfGspaceLoop
	(	int iGsq = iG.length_squared();
		int iGz = iG[2];
		if(iGsq == iGz * iGz) // => iG || z (Cauchy-Schwarz)
		{	double GzH = (2*M_PI) * fabs(zH * iGz);
			double GzSigma = zSigma * sqrt(GGT.metric_length_squared(iG));
			w[i] = prefac * bessel_jl(0,GzH) //norm and width
				* cis(-(2.*M_PI)*iGz*z0) //translate
				* exp(-0.5*GzSigma*GzSigma); //Gauss smoothing
		}
		else w[i] = 0.;
	)
}

//Save slab-weight matrix elements in Wannier basis, if requested:
void WannierMinimizer::saveMLWF_W(int iSpin, const matrix& phase)
{	assert(wannier.zH);
	//Construct slab weight function:
	ScalarFieldArray w(1);
	{	ScalarFieldTilde wTilde; nullToZero(wTilde, e.gInfo);
		threadLaunch(slabWeight_thread, e.gInfo.nG, e.gInfo.S, e.gInfo.GGT,
			wTilde->data(), wannier.z0, wannier.zH, wannier.zSigma);
		w[0] = I(wTilde);
	}
	saveMLWF(iSpin, phase, w, "mlwfW", false);
}

//Thread function for setting z position operator:
inline void zWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, double Lz, double z0, double* z)
{	double invSz = 1./S[2];
	THREAD_rLoop
	(	z[i] = invSz*iv[2] - z0;
		z[i] -= floor(0.5 + z[i]); //wrap to [-0.5,0.5)
		if(fabs(z[i]) > 0.5-symmThreshold) z[i] = 0.; //Nyquist component
		z[i] *= Lz; //convert to Cartesian bohrs
	)
}

void readDensityArray(ScalarFieldArray& var, string varName, string fnamePattern, const Everything* e); //declared in ElecVars.cpp

//Save z position matrix elements in Wannier basis, if requested:
void WannierMinimizer::saveMLWF_Z(int iSpin, const matrix& phase)
{	assert(wannier.zVfilename.length());
	ScalarFieldArray z; 
	if(wannier.zVfilename == "Ramp")
	{	//Construct ramp along third lattice direction:
		assert(e.coulombParams.isTruncated()[2]);
		nullToZero(z, e.gInfo, 1);
		double Lz = e.gInfo.R.column(2).length();
		double z0 = e.coulombParams.embedCenter[2]; //in lattice coordinates
		threadLaunch(zWeight_thread, e.gInfo.nr, e.gInfo.S, Lz, z0, z[0]->data());
		z[0] = I(gaussConvolve(J(z[0]), 1.)); //smooth transition at boundary
	}
	else
	{	//Construct from difference in potentials (account for local field):
		assert(wannier.zFieldMag); //must have calculation with a different field read in
		ScalarFieldArray V0, V1;
		readDensityArray(V0, "Vscloc", wannier.initFilename, &e); //current V
		readDensityArray(V1, "Vscloc", wannier.zVfilename, &e); //V from changed Efield
		z = (1./(wannier.zFieldMag*e.gInfo.dV)) * (V1-V0);
	}
	saveMLWF(iSpin, phase, z, "mlwfZ", true);
}

//Helper function for scalar field matrix elements (eg. slab and z)
void WannierMinimizer::saveMLWF(int iSpin, const matrix& phase, const ScalarFieldArray& w, string varName, bool suppressUnbound)
{	//Prepare scalar field for direct and augmented contributions:
	ScalarFieldArray JdagOJw;
	for(const ScalarField& wi: w)
		JdagOJw.push_back(JdagOJ(wi));
	e.iInfo.augmentDensityGridGrad(JdagOJw);
	//Compute matrix elements of w between Bloch states at full k mesh:
	matrix wWannierTilde = zeroes(nCenters*nCenters, phase.nRows());
	int iqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	//Direct contributions from bands:
		const ColumnBundle& Ck = getWfns(kMesh[i].point, iSpin);
		matrix wSub = Ck ^ Idag_DiagV_I(Ck, JdagOJw);
		//Ultrasoft augmentation:
		for(const auto& sp: e.iInfo.species)
			if(sp->isUltrasoft())
			{	matrix VdagCk = (*sp->getV(Ck)) ^ Ck, wVdagCk;
				sp->augmentDensitySphericalGrad(*(Ck.qnum), VdagCk, wVdagCk);
				wSub += dagger(VdagCk) * wVdagCk;
			}
		//Suppress matrix elements of unbound states for numerical stability:
		if(suppressUnbound and e.coulombParams.geometry!=CoulombParams::Periodic)
		{	diagMatrix weightE = e.eVars.Hsub_eigs[kMesh[i].point.iReduced];
			double sigmaE = 0.02; //selected to transition weights over ~ 1 eV below vacuum level
			for(double& wE: weightE) wE = 0.5*erfc(wE/sigmaE+1); //convert energy to weight -> 0 for unbound states
			wSub = weightE * wSub * weightE;
		}
		wSub = dagger(kMesh[i].U) * wSub * kMesh[i].U; //apply MLWF-optimized rotations
		callPref(eblas_copy)(wWannierTilde.dataPref()+wWannierTilde.index(0,iqMine), wSub.dataPref(), wSub.nData());
		iqMine++;
	}
	//Fourier transform to Wannier space and save
	dumpWannierized(wWannierTilde, phase, varName, realPartOnly, iSpin);
}


//Electron-electron linewidths:
void WannierMinimizer::saveMLWF_ImSigma_ee(int iSpin, const matrix& phase)
{	string fname = wannier.getFilename(Wannier::FilenameInit, "ImSigma_ee");
	if(fileSize(fname.c_str()) >= 0)
	{	//Read Bloch version:
		logPrintf("Reading '%s' ... ", fname.c_str()); logFlush();
		std::vector<diagMatrix> ImSigma_ee;
		e.eInfo.read(ImSigma_ee, fname.c_str());
		logPrintf("done.\n"); logFlush();
		//Fill in linewidths for states out of range (to ensure smoothness in Wannier interpolation)
		//--- determine weights in fit:
		std::vector<diagMatrix> fitWeight(e.eInfo.nStates);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	fitWeight[q].reserve(nBands);
			for(double IS: ImSigma_ee[q]) //first pass: whether ImSigma is valid or not
				fitWeight[q].push_back(std::isnan(IS) ? 0. : 1.);
			for(int b=0; b<nBands; b++) if(fitWeight[q][b])
			{	double minEdist = 100.; //keep finite so that no error below if all ImSigma's already defined
				for(int b2=0; b2<nBands; b2++)
					if(!fitWeight[q][b2])
						minEdist = std::min(minEdist, fabs(e.eVars.Hsub_eigs[q][b]-e.eVars.Hsub_eigs[q][b2]));
				fitWeight[q][b] = 1./hypot(minEdist, 0.1); //cap the weights at some reasonable amount
			}
		}
		//--- construct matrices for polynomial fit:
		int order = 3; //quadratic
		int nRows = e.eInfo.nStates * nBands;
		matrix Lhs = zeroes(nRows, order);
		matrix rhs = zeroes(nRows, 1);
		diagMatrix w(nRows, 0.);
		int row = e.eInfo.qStart * nBands;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			for(int b=0; b<nBands; b++)
			{	double E = e.eVars.Hsub_eigs[q][b], Epow = 1.;
				for(int n=0; n<order; n++)
				{	Lhs.set(row,n, Epow);
					Epow *= E;
				}
				rhs.set(row,0, fitWeight[q][b] ? ImSigma_ee[q][b] : 0.); //Handle NANs correctly
				w[row] = fitWeight[q][b];
				row++;
			}
		mpiWorld->allReduceData(Lhs, MPIUtil::ReduceSum);
		mpiWorld->allReduceData(rhs, MPIUtil::ReduceSum);
		mpiWorld->allReduceData(w, MPIUtil::ReduceSum);
		//--- weighted least squares polynomial fit
		matrix rhsFit = Lhs * (inv(dagger(Lhs)*w*Lhs) * (dagger(Lhs)*w*rhs));
		//--- fill in missing values
		row = e.eInfo.qStart * nBands;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			for(int b=0; b<nBands; b++)
			{	if(std::isnan(ImSigma_ee[q][b]))
					ImSigma_ee[q][b] =  rhsFit(row,0).real();
				row++;
			}
		//Wannierize:
		matrix ImSigma_eeWannierTilde = zeroes(nCenters*nCenters, phase.nRows());
		int iqMine = 0;
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	matrix ImSigma_eeSub = dagger(kMesh[i].U) * ImSigma_ee[kMesh[i].point.iReduced + iSpin*qCount] * kMesh[i].U;
			callPref(eblas_copy)(ImSigma_eeWannierTilde.dataPref()+ImSigma_eeWannierTilde.index(0,iqMine), ImSigma_eeSub.dataPref(), ImSigma_eeSub.nData());
			iqMine++;
		}
		dumpWannierized(ImSigma_eeWannierTilde, phase, "mlwfImSigma_ee", realPartOnly, iSpin);
	}
}
