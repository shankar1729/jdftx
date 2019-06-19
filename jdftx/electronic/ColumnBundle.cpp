/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/matrix.h>
#include <core/vector3.h>
#include <core/Random.h>
#include <core/BlasExtra.h>
#include <core/ScalarFieldIO.h>
#include <fftw3.h>

// Called by other constructors to do the work
void ColumnBundle::init(int nc, size_t len, const Basis *b, const QuantumNumber* q, bool onGpu)
{
	ncols = nc;
	col_length = len;
	basis = b;
	qnum = q;

	if(nCols() == 0) { memFree(); return; } //must be default constructor or assignment to empty ColumnBundle
	assert(colLength() != 0);
	memInit("ColumnBundle", nCols()*colLength(), onGpu); //in base class ManagedMemory
}

void ColumnBundle::free()
{	ncols = 0;
	col_length = 0;
	basis = 0;
	qnum = 0;
	memFree();
}

// (Default) Constructor
ColumnBundle::ColumnBundle(int nc, size_t len, const Basis *b, const QuantumNumber* q, bool onGpu)
{	init(nc, len, b, q, onGpu);
}
// Copy constructor
ColumnBundle::ColumnBundle(const ColumnBundle &Y)
{	init(Y.nCols(), Y.colLength(), Y.basis, Y.qnum, Y.isOnGpu()); //initialize size and storage
	if(nData()) memcpy((ManagedMemory&)*this, (const ManagedMemory&)Y); //copy data
}
// Move constructor
ColumnBundle::ColumnBundle(ColumnBundle&& Y)
{	std::swap(ncols, Y.ncols);
	std::swap(col_length, Y.col_length);
	std::swap(qnum, Y.qnum);
	std::swap(basis, Y.basis);
	memMove((ManagedMemory&&)Y); //cannibalize Y's data
}
// Create a similar object, without copying data
ColumnBundle ColumnBundle::similar(int ncOverride) const
{	return ColumnBundle(ncOverride<0 ? ncols : ncOverride, col_length, basis, qnum, isOnGpu());
}


// Copy-assignment
ColumnBundle& ColumnBundle::operator=(const ColumnBundle &Y)
{	init(Y.nCols(), Y.colLength(), Y.basis, Y.qnum, Y.isOnGpu()); //initialize size and storage
	if(nData()) memcpy((ManagedMemory&)*this, (const ManagedMemory&)Y); //copy data
	return *this;
}
// Move-assignment
ColumnBundle& ColumnBundle::operator=(ColumnBundle&& Y)
{	std::swap(ncols, Y.ncols);
	std::swap(col_length, Y.col_length);
	std::swap(qnum, Y.qnum);
	std::swap(basis, Y.basis);
	memMove((ManagedMemory&&)Y); //cannibalize Y's data
	return *this;
}
ColumnBundle clone(const ColumnBundle& Y)
{	return Y;
}

void randomize(ColumnBundle& x)
{	x.randomize(0, x.nCols());
}

double dot(const ColumnBundle& x, const ColumnBundle& y)
{	complex result = dotc(x, y)*2.0;
	return result.real();
}

ColumnBundle ColumnBundle::getSub(int colStart, int colStop) const
{	assert(colStart>=0);
	assert(colStop<=nCols());
	int nColsSub = colStop - colStart;
	assert(nColsSub>0);
	ColumnBundle ret = this->similar(nColsSub);
	callPref(eblas_copy)(ret.dataPref(), dataPref()+colStart*colLength(), nColsSub*colLength());
	return ret;
}

void ColumnBundle::setSub(int colStart, const ColumnBundle& Y)
{	assert(colStart>=0);
	assert(colStart<nCols());
	assert(colLength()==Y.colLength());
	int nColsSub = std::min(Y.nCols(), nCols()-colStart);
	callPref(eblas_copy)(dataPref()+colStart*colLength(), Y.dataPref(), nColsSub*colLength());
}

#define CHECK_COLUMN_INDEX \
	assert(i>=0 && i<nCols()); \
	assert(s>=0 && s<spinorLength());

complexScalarFieldTilde ColumnBundle::getColumn(int i, int s) const
{	const GridInfo& gInfo = *(basis->gInfo);
	CHECK_COLUMN_INDEX
	complexScalarFieldTilde full; nullToZero(full, gInfo); //initialize a full G-space vector to zero
	//scatter from the i'th column to the full vector:
	callPref(eblas_scatter_zdaxpy)(basis->nbasis, 1., basis->index.dataPref(), dataPref()+index(i,s*basis->nbasis), full->dataPref());
	return full;
}

void ColumnBundle::setColumn(int i, int s, const complexScalarFieldTilde& full)
{	//Zero the i'th column:
	callPref(eblas_zero)(basis->nbasis, dataPref()+index(i,s*basis->nbasis));
	//Gather-accumulate from the full vector into the i'th column
	accumColumn(i,s, full);
}

void ColumnBundle::accumColumn(int i, int s, const complexScalarFieldTilde& full)
{	assert(full);
	CHECK_COLUMN_INDEX
	//Gather-accumulate from the full vector into the i'th column
	callPref(eblas_gather_zdaxpy)(basis->nbasis, 1., basis->index.dataPref(), full->dataPref(), dataPref()+index(i,s*basis->nbasis));
}
#undef CHECK_COLUMN_INDEX


// Allocate an array of ColumnBundles
void init(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const Basis* basis, const ElecInfo* eInfo)
{	Y.resize(nbundles);
	if(ncols && basis && eInfo)
	{	assert(nbundles >= eInfo->qStop);
		for(int q=eInfo->qStart; q<eInfo->qStop; q++)
			Y[q].init(ncols, basis[q].nbasis * eInfo->spinorLength(), basis+q, &eInfo->qnums[q], isGpuEnabled());
	}
}


// Randomize with a high frequency cutoff of 0.75 hartrees
void ColumnBundle::randomize(int colStart, int colStop)
{	static StopWatch watch("ColumnBundle::randomize"); watch.start();
	assert(basis->nbasis==colLength() || 2*basis->nbasis==colLength());
	complex* thisData = data(); //currently only on cpu
	int nSpinor = colLength()/basis->nbasis;
	size_t j=0; //basis index
	for(const vector3<int>& iG: basis->iGarr)
	{	vector3<> kplusG = iG + qnum->k;
		double KE = 0.5*dot(kplusG, basis->gInfo->GGT*kplusG);
		double t = KE/0.75;
		double sigma = 1.0/((1.0+t*t*t*t*t*t) * basis->gInfo->detR);
		for(int s=0; s<nSpinor; s++)
			for(int i=colStart; i<colStop; i++)
				thisData[index(i,j+s*basis->nbasis)] = Random::normalComplex(sigma);
		j++;
	}
	watch.stop();
}
void randomize(std::vector<ColumnBundle>& Y, const ElecInfo& eInfo)
{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		if(Y[q]) Y[q].randomize(0, Y[q].nCols());
}

//--------- Read/write an array of ColumnBundles from/to a file --------------

void ElecInfo::write(const std::vector<ColumnBundle>& Y, const char* fname) const
{
#if MPI_SAFE_WRITE
	//Safe mode / write from head:
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname, "w");
		if(!fp) die_alone("Error opening file '%s' for writing.\n", fname);
		for(int q=0; q<nStates; q++)
		{	if(!isMine(q))
			{	size_t nData = 0;
				mpiWorld->recv(nData, whose(q), q);
				ManagedArray<complex> buf; buf.init(nData);
				mpiWorld->recvData(buf, whose(q), q);
				buf.write(fp);
			}
			else Y[q].write(fp);
		}
		fclose(fp);
	}
	else
		for(int q=qStart; q<qStop; q++)
		{	mpiWorld->send(Y[q].nData(), 0, q);
			mpiWorld->sendData(Y[q], 0, q);
		}
#else
	//Compute output length from each process:
	std::vector<long> nBytes(mpiWorld->nProcesses(), 0); //total bytes to be written on each process
	for(int q=qStart; q<qStop; q++)
		nBytes[mpiWorld->iProcess()] += Y[q].nData()*sizeof(complex);
	//Sync nBytes across processes:
	if(mpiWorld->nProcesses()>1)
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
			mpiWorld->bcast(nBytes[iSrc], iSrc);
	//Compute offset of current process, and expected file length:
	long offset=0, fsize=0;
	for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
	{	if(iSrc<mpiWorld->iProcess()) offset += nBytes[iSrc];
		fsize += nBytes[iSrc];
	}
	//Write to file:
	MPIUtil::File fp; mpiWorld->fopenWrite(fp, fname);
	mpiWorld->fseek(fp, offset, SEEK_SET);
	for(int q=qStart; q<qStop; q++)
		mpiWorld->fwriteData(Y[q], fp);
	mpiWorld->fclose(fp);
#endif
}


ElecInfo::ColumnBundleReadConversion::ColumnBundleReadConversion()
: realSpace(false), nBandsOld(0), Ecut(0), EcutOld(0)
{
}

void ElecInfo::read(std::vector<ColumnBundle>& Y, const char *fname, const ColumnBundleReadConversion* conversion) const
{	if(conversion && conversion->realSpace)
	{	if(qStop==qStart) return; //no k-point on this process
		const GridInfo* gInfoWfns = Y[qStart].basis->gInfo;
		//Create a custom gInfo if necessary:
		GridInfo gInfoCustom;
		gInfoCustom.R = gInfoWfns->R;
		gInfoCustom.S = conversion->S_old;
		for(int k=0; k<3; k++) if(!gInfoCustom.S[k]) gInfoCustom.S[k] = gInfoWfns->S[k];
		bool needCustom = !(gInfoCustom.S == gInfoWfns->S);
		if(needCustom) { logSuspend(); gInfoCustom.initialize(); logResume(); }
		const GridInfo& gInfo = needCustom ? gInfoCustom : *gInfoWfns;
		//Read one column at a time:
		complexScalarField Icol; nullToZero(Icol, gInfo);
		for(int q=qStart; q<qStop; q++)
		{	int nCols = Y[q].nCols();
			int nSpinor = Y[q].spinorLength();
			if(conversion->nBandsOld) nCols = std::min(nCols, conversion->nBandsOld);
			for(int b=0; b<nCols; b++) for(int s=0; s<nSpinor; s++)
			{	char fname_qb[1024]; sprintf(fname_qb, fname, q, b*nSpinor+s);
				loadRawBinary(Icol, fname_qb);
				if(needCustom) Y[q].setColumn(b,s, changeGrid(J(Icol), *gInfoWfns));
				else Y[q].setColumn(b,s, J(Icol));
			}
		}
	}
	else
	{	//Check if a conversion is actually needed:
		std::vector<ColumnBundle> Ytmp(qStop);
		std::vector<Basis> basisTmp(qStop);
		std::vector<long> nBytes(mpiWorld->nProcesses(), 0); //total bytes to be read on each process
		for(int q=qStart; q<qStop; q++)
		{	bool needTmp = false, customBasis = false;
			int nCols = Y[q].nCols();
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
					basisTmp[q].setup(*(Y[q].basis->gInfo), *(Y[q].basis->iInfo), EcutOld, Y[q].qnum->k);
					logResume();
				}
			}
			const Basis* basis = customBasis ? &basisTmp[q] : Y[q].basis;
			int nSpinor = Y[q].spinorLength();
			if(needTmp) Ytmp[q].init(nCols, basis->nbasis*nSpinor, basis, Y[q].qnum);
			nBytes[mpiWorld->iProcess()] += nCols * basis->nbasis*nSpinor * sizeof(complex);
		}
		//Sync nBytes:
		if(mpiWorld->nProcesses()>1)
			for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
				mpiWorld->bcast(nBytes[iSrc], iSrc);
		//Compute offset of current process, and expected file length:
		long offset=0, fsize=0;
		for(int iSrc=0; iSrc<mpiWorld->nProcesses(); iSrc++)
		{	if(iSrc<mpiWorld->iProcess()) offset += nBytes[iSrc];
			fsize += nBytes[iSrc];
		}
		//Read data into Ytmp or Y as appropriate, and convert if necessary:
		MPIUtil::File fp; mpiWorld->fopenRead(fp, fname, fsize,
			(e->vibrations and qnums.size()>1)
			? "Hint: Vibrations requires wavefunctions without symmetries:\n"
				"either don't read in state, or consider using phonon instead.\n"
			: "Hint: Did you specify the correct nBandsOld, EcutOld and kdepOld?\n");
		mpiWorld->fseek(fp, offset, SEEK_SET);
		for(int q=qStart; q<qStop; q++)
		{	ColumnBundle& Ycur = Ytmp[q] ? Ytmp[q] : Y[q];
			mpiWorld->freadData(Ycur, fp);
			if(Ytmp[q]) //apply conversions:
			{	if(Ytmp[q].basis!=Y[q].basis)
				{	int nSpinor = Y[q].spinorLength();
					for(int b=0; b<std::min(Y[q].nCols(), Ytmp[q].nCols()); b++)
						for(int s=0; s<nSpinor; s++)
							Y[q].setColumn(b,s, Ytmp[q].getColumn(b,s)); //convert using the full G-space as an intermediate
				}
				else
				{	if(Ytmp[q].nCols()<Y[q].nCols()) Y[q].setSub(0, Ytmp[q]);
					else Y[q] = Ytmp[q].getSub(0, Y[q].nCols());
				}
				Ytmp[q].free();
			}
		}
		mpiWorld->fclose(fp);
	}
}
