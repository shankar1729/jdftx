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
#include <electronic/matrix.h>
#include <core/vector3.h>
#include <core/Random.h>
#include <core/BlasExtra.h>
#include <fftw3.h>

// Called by other constructors to do the work
void ColumnBundle::init(int nc, size_t len, const Basis *b, const QuantumNumber* q, bool onGpu)
{
	ncols = nc;
	col_length = len;
	basis = b;
	qnum = q;

	if(nCols() == 0) return; //must be default constructor
	assert(colLength() != 0);
	memInit(nCols()*colLength(), onGpu); //in base class ManagedMemory
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
{	reader = 0;
	init(nc, len, b, q, onGpu);
}
// Copy constructor
ColumnBundle::ColumnBundle(const ColumnBundle &Y)
{	init(Y.nCols(), Y.colLength(), Y.basis, Y.qnum, Y.isOnGpu()); //initialize size and storage
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)Y); //copy data
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
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)Y); //copy data
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


complexDataGptr ColumnBundle::getColumn(int i) const
{	const GridInfo& gInfo = *(basis->gInfo);
	complexDataGptr full; nullToZero(full, gInfo); //initialize a full G-space vector to zero
	//scatter from the i'th column to the full vector:
	callPref(eblas_scatter_zdaxpy)(basis->nbasis, 1, basis->indexPref, dataPref()+index(i,0), full->dataPref());
	return full;
}

void ColumnBundle::setColumn(int i, const complexDataGptr& full)
{	//Zero the i'th column:
	callPref(eblas_zero)(basis->nbasis, dataPref()+index(i,0));
	//Gather-accumulate from the full vector into the i'th column
	accumColumn(i, full);
}

void ColumnBundle::accumColumn(int i, const complexDataGptr& full)
{	//Gather-accumulate from the full vector into the i'th column
	callPref(eblas_gather_zdaxpy)(basis->nbasis, 1.0, basis->indexPref, full->dataPref(), dataPref()+index(i,0));
}



// Allocate an array of ColumnBundles
void init(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const Basis* basis, const QuantumNumber* qnum)
{	Y.resize(nbundles);
	if(ncols && basis && qnum)
		for(int i=0; i<nbundles; i++)
			Y[i].init(ncols, basis[i].nbasis, basis+i, qnum+i, isGpuEnabled());
}


// Randomize with a high frequency cutoff of 0.75 hartrees
void ColumnBundle::randomize(int colStart, int colStop)
{	assert(basis->nbasis == colLength());
	complex* thisData = data(); //currently only on cpu
	for(size_t j=0; j < basis->nbasis; j++)
	{	vector3<> kplusG = basis->iGarr[j] + qnum->k;
		double KE = 0.5*dot(kplusG, basis->gInfo->GGT*kplusG);
		double t = KE/0.75;
		double sigma = 1.0/(1.0+t*t*t*t*t*t);
		for(int i=colStart; i < colStop; i++)
			thisData[index(i,j)] = Random::normalComplex(sigma);
	}
}
void randomize(std::vector<ColumnBundle>& Y)
{	for(ColumnBundle& y: Y) y.randomize(0, y.nCols());
}

void write(const std::vector<ColumnBundle>& Y, const char* fname)
{	FILE *fp = fopen(fname, "w");
	if(!fp) die("Error opening %s for writing.\n", fname);
	for(const ColumnBundle& y: Y) y.write(fp);
	fclose(fp);
}


//-------- Column bundle read() with grid, band change and real-space options ----------

struct ColumnBundleReader
{	int* indexMap; //mapping from old to new indices: i.e. newIndex = indexMap[oldIndex];
	size_t lengthOld, length; //old and new colum lengths
	int nBandsOld, nBands; //old and new band number

	ColumnBundleReader(const Everything& e, const ColumnBundle& cb):indexMap(0)
	{	nBandsOld = e.eVars.nBandsOld; nBands = e.eInfo.nBands;
		if(e.eVars.EcutOld != e.cntrl.Ecut
			|| e.eVars.kdepOld != e.cntrl.basisKdep)
		{	//Need to interpolate
			const Basis& b = *cb.basis;
			length = b.nbasis;
			int ibox[3];
			for(int k=0; k<3; k++)
				ibox[k] = 1 + int(e.gInfo.R.column(k).length()
					* sqrt(2.0*std::max(e.eVars.EcutOld, e.cntrl.Ecut)) / (2*M_PI));
			register matrix3<> GGT = e.gInfo.GGT;
			//Calculate input length
			vector3<> kOld = (e.eVars.kdepOld==BasisKpointDep ? cb.qnum->k : vector3<>(0,0,0));
			lengthOld = 0;
			for (int i0 = -ibox[0]; i0 <= ibox[0]; i0++)
				for (int i1 = -ibox[1]; i1 <= ibox[1]; i1++)
					for (int i2 = -ibox[2]; i2 <= ibox[2]; i2++)
					{
						register vector3<> fOld(i0,i1,i2); fOld += kOld;
						if(0.5*dot(fOld, GGT*fOld) <= e.eVars.EcutOld)
							lengthOld++;
					}
			//Create index map:
			indexMap = new int[lengthOld];
			register int indexOld=0, index=0;
			vector3<> k = (e.cntrl.basisKdep==BasisKpointDep ? cb.qnum->k : vector3<>(0,0,0));
			for (int i0 = -ibox[0]; i0 <= ibox[0]; i0++)
				for (int i1 = -ibox[1]; i1 <= ibox[1]; i1++)
					for (int i2 = -ibox[2]; i2 <= ibox[2]; i2++)
					{
						register vector3<> f(i0,i1,i2); f+=k; register bool in = 0.5*dot(f,GGT*f) <= e.cntrl.Ecut;
						register vector3<> fOld(i0,i1,i2); fOld += kOld; register bool inOld = 0.5*dot(fOld,GGT*fOld) <= e.eVars.EcutOld;
						if(inOld)
						{	indexMap[indexOld] = (in ? index : -1); //negative index points discarded during transfer
							indexOld++;
						}
						if(in) index++;
					}
		}
		else
			length = lengthOld = cb.basis->nbasis;
	}
	~ColumnBundleReader()
	{	if(indexMap) delete[] indexMap;
	}

	void read(FILE* fp, ColumnBundle& cb)
	{	complex* buf = new complex[lengthOld];
		cb.zero();
		for(int b=0; b<nBandsOld; b++)
		{	complex* cb_b_data = cb.data() + cb.index(b,0);
			fread(buf,sizeof(complex),lengthOld,fp);
			if(b<nBands)
			{
				if(indexMap)
				{	//Need conversion:
					for(size_t i=0; i<lengthOld; i++)
						if(indexMap[i]>=0)
							cb_b_data[indexMap[i]] = buf[i];
				}
				else
				{	//No conversion:
					memcpy(cb_b_data, buf, sizeof(complex)*length);
				}
			}
		}
		delete[] buf;
	}
};

void ColumnBundle::read(FILE *fp)
{	reader->read(fp, *this);
}

size_t ColumnBundle::createReader(const Everything& e)
{	reader = new ColumnBundleReader(e, *this);
	return reader->nBandsOld * reader->lengthOld * sizeof(complex);
}
void ColumnBundle::destroyReader()
{	delete reader;
}

void readRealSpace(std::vector<ColumnBundle>& Y, const char *fnamePattern, const Everything& e)
{	const ElecVars& eVars = e.eVars;
	const ElecInfo& eInfo = e.eInfo;
	int nBandsOld = eVars.nBandsOld, nBands = eInfo.nBands;
	int Nx = eVars.NxOld;
	int Ny = eVars.NyOld;
	int Nz = eVars.NzOld;
	size_t lengthOld = Nx * Ny * Nz;
	complex* buf = (complex*)fftw_malloc(lengthOld*sizeof(complex));
	fftw_plan planJ = fftw_plan_dft_3d(Nx,Ny,Nz, (fftw_complex*)buf,(fftw_complex*)buf, FFTW_FORWARD, FFTW_MEASURE);
	double scaleFac = sqrt(e.gInfo.detR)/lengthOld;

	for(unsigned q=0; q<Y.size(); q++)
	{	const Basis& b = *Y[q].basis;
		size_t length = b.nbasis;
		int* indexMap = new int[length];

		//Initialize the index map:
		register size_t index=0;
		const matrix3<>& GGT = e.gInfo.GGT;
		vector3<> k = (e.cntrl.basisKdep==BasisKpointDep ? Y[q].qnum->k : vector3<>(0,0,0));
		int ibox[3]; for(int i=0; i<3; i++) ibox[i] = (b.gInfo->S[i]-2)/4+1;
		for (int i0 = -ibox[0]; i0 <= ibox[0]; i0++)
			for (int i1 = -ibox[1]; i1 <= ibox[1]; i1++)
				for (int i2 = -ibox[2]; i2 <= ibox[2]; i2++)
				{
					bool inOld = true; register int indexOld=0;
					if(abs(i0)<Nx/2) indexOld = indexOld*Nx + (i0>=0 ? i0 : (Nx+i0)); else inOld=false;
					if(abs(i1)<Ny/2) indexOld = indexOld*Ny + (i1>=0 ? i1 : (Ny+i1)); else inOld=false;
					if(abs(i2)<Nz/2) indexOld = indexOld*Nz + (i2>=0 ? i2 : (Nz+i2)); else inOld=false;

					register vector3<> f(i0,i1,i2); f += k;
					if(0.5*dot(f,GGT*f) <= e.cntrl.Ecut)
						indexMap[index++] = (inOld ? indexOld : -1); //negative index points discarded during transfer
				}
		assert(index==length);

		for(int b=0; b<std::min(nBands,nBandsOld); b++)
		{	char fname[1024]; sprintf(fname, fnamePattern, q, b);
			FILE* fp = fopen(fname, "rb");
			size_t nRead = fread(buf, sizeof(complex), lengthOld, fp);
			if(nRead < lengthOld) die("File '%s' ended too soon: %lu records of expected %lu.\n", fname, nRead, lengthOld);

			fftw_execute(planJ); //transforms buf in-place
			complex* Yq_b_data = Y[q].data() + Y[q].index(b,0);
			for(size_t i=0; i<length; i++)
				 Yq_b_data[i] = (indexMap[i]>=0 ? scaleFac*buf[indexMap[i]] : 0.0);
		}
		delete[] indexMap;
	}
	delete[] buf;
	fftw_destroy_plan(planJ);
}

//Get the size of a file:
#include <sys/stat.h>
off_t fsize(const char *filename)
{	struct stat st;
	if(stat(filename, &st) == 0) return st.st_size;
    return -1;
}


// Read/write an array of ColumnBundles from/to a file
void read(std::vector<ColumnBundle>& Y, const char *fname, const Everything& e)
{
	if(e.eVars.readWfnsRealspace) readRealSpace(Y, fname, e);
	else
	{	//Fourier space read:
		off_t fLen = fsize(fname);
		if(fLen<0) die("Error accessing '%s'.\n", fname);

		//make sure file size exactly matches expectation
		off_t expectedLen=0;
		for(ColumnBundle& y: Y) expectedLen += y.createReader(e);
		if(expectedLen!=fLen)
		{	die("Length of '%s' was %ld instead of the expected %ld bytes.\n"
				"Hint: Did you specify the correct nBandsOld, EcutOld and kdepOld?\n",
				fname, fLen, expectedLen);
		}

		FILE *fp = fopen(fname,"rb");
		for(ColumnBundle& y: Y)
		{	y.read(fp);
			y.destroyReader();
		}
		fclose(fp);
	}
}

//------------------------ Arithmetic --------------------

ColumnBundle& operator+=(ColumnBundle& Y, const ColumnBundle &X) { if(Y) axpy(1.0, X, Y); else Y=X; return Y; }
ColumnBundle& operator-=(ColumnBundle& Y, const ColumnBundle &X) { if(Y) axpy(-1.0, X, Y); else Y=-X; return Y; }
ColumnBundle operator+(const ColumnBundle &Y1, const ColumnBundle &Y2) { ColumnBundle Ysum(Y1); Ysum += Y2; return Ysum; }
ColumnBundle operator-(const ColumnBundle &Y1,const ColumnBundle &Y2) { ColumnBundle Ydiff(Y1); Ydiff -= Y2; return Ydiff; }

ColumnBundle& operator*=(ColumnBundle& X, double s) { scale(s, X); return X; }
scaled<ColumnBundle> operator*(double s, const ColumnBundle &Y) { return scaled<ColumnBundle>(Y, s); }
scaled<ColumnBundle> operator*(const ColumnBundle &Y, double s) { return scaled<ColumnBundle>(Y, s); }
scaled<ColumnBundle> operator-(const ColumnBundle &Y) { return scaled<ColumnBundle>(Y, -1); }
ColumnBundle& operator*=(ColumnBundle& X, complex s) { scale(s, X); return X; }
ColumnBundle operator*(complex s, const ColumnBundle &Y) { ColumnBundle sY(Y); sY *= s; return sY; }
ColumnBundle operator*(const ColumnBundle &Y, complex s) { ColumnBundle sY(Y); sY *= s; return sY; }

ColumnBundle operator*(const scaled<ColumnBundle> &sY, const matrixScaledTransOp &Mst)
{	static StopWatch watch("Y*M");
	watch.start();
	const ColumnBundle& Y = sY.data;
	assert(Y.nCols()==Mst.nRows());
	const matrix& M = Mst.mat;
	double scaleFac = sY.scale * Mst.scale;
	ColumnBundle YM = Y.similar(Mst.nCols());
	callPref(eblas_zgemm)(CblasNoTrans, Mst.op, YM.colLength(), YM.nCols(), Y.nCols(),
		scaleFac, Y.dataPref(), Y.colLength(), M.dataPref(), M.nRows(),
		0.0, YM.dataPref(), YM.colLength());
	watch.stop();
	return YM;
}
ColumnBundle operator*(const scaled<ColumnBundle> &sY, const diagMatrix& d)
{	const ColumnBundle& Y = sY.data;
	assert(Y.nCols()==d.nRows());
	ColumnBundle Yd = Y; complex* YdData = Yd.dataPref();
	for(int i=0; i<d.nCols(); i++)
		callPref(eblas_zscal)(Yd.colLength(), sY.scale*d[i], YdData+Yd.index(i,0), 1);
	return Yd;
}
matrix operator^(const scaled<ColumnBundle> &sY1, const scaled<ColumnBundle> &sY2)
{	static StopWatch watch("Y1^Y2");
	watch.start();
	const ColumnBundle& Y1 = sY1.data;
	const ColumnBundle& Y2 = sY2.data;
	double scaleFac = sY1.scale * sY2.scale;
	assert(Y1.colLength() == Y2.colLength());
	matrix Y1dY2(Y1.nCols(), Y2.nCols(), isGpuEnabled());
	callPref(eblas_zgemm)(CblasConjTrans, CblasNoTrans, Y1.nCols(), Y2.nCols(), Y1.colLength(),
		scaleFac, Y1.dataPref(), Y1.colLength(), Y2.dataPref(), Y2.colLength(),
		0.0, Y1dY2.dataPref(), Y1dY2.nRows());
	watch.stop();
	return Y1dY2;
}
