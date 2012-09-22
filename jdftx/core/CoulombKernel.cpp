/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <core/CoulombKernel.h>
#include <core/Coulomb_internal.h>
#include <core/LatticeUtils.h>
#include <core/LoopMacros.h>
#include <core/Util.h>
#include <core/Thread.h>

const double CoulombKernelDesc::nSigmasPerWidth = 10.;


bool isOrthogonal(const vector3<>& a, const vector3<>& b)
{	return pow(dot(a, b), 2) < symmThresholdSq * a.length_squared() * b.length_squared();
}

CoulombKernelDesc::CoulombKernelDesc(const matrix3<> R, const vector3<int> S,
	const vector3<bool> isTruncated, const vector3<> sigmaBorder, double omega)
: R(R), S(S), isTruncated(isTruncated), sigmaBorder(sigmaBorder), omega(omega)
{
	//For each pair of directions:
	std::vector<std::pair<int,int>> brokenSym; //pairs of directions for which symmetries are broken by truncation geometry
	for(int iDir=0; iDir<3; iDir++)
	{	int jDir = (iDir+1)%3;
		//make sure non-orthogonal ones have same sigmaBorder
		if(!isOrthogonal(R.column(iDir), R.column(jDir)))
			assert(sigmaBorder[iDir] == sigmaBorder[jDir]);
		if( (isTruncated[iDir]!=isTruncated[jDir])
			|| (sigmaBorder[iDir]!=sigmaBorder[jDir]) )
			brokenSym.push_back(std::make_pair(iDir, jDir));
	}
	
	//Compute symmetries for kernel I/O compression:
	std::vector<matrix3<int>> symLattice = getSymmetries(R); //symmetries of the Bravais lattice
	for(const matrix3<int>& m: symLattice)
	{	//Check if matrix connects directions whose symmetry is broken by truncation geometry
		bool broken = false;
		for(std::pair<int,int> d: brokenSym)
			if(m(d.first,d.second) || m(d.second,d.first))
			{	broken = true;
				break;
			}
		if(!broken) sym.push_back(m);
	}
}

//! Loop over the symmetric images of a reciprocal lattice vector iG
#define LOOP_symImagesG(code) \
	for(const matrix3<int>& m: sym) \
	{	vector3<int> jG = iG * m; /* image point (note that rotation is by transpose of m) */ \
		/* Check boundaries and wrap to positive domain: */ \
		if(abs(2*jG[0]) > S[0]) continue; if(jG[0]<0) jG[0] += S[0]; \
		if(abs(2*jG[1]) > S[1]) continue; if(jG[1]<0) jG[1] += S[1]; \
		if(abs(2*jG[2]) > S[2]) continue; if(jG[2]<0) continue; \
		int j = jG[2] + size2*(jG[1] + S[1]*jG[0]); \
		code \
	}

static const size_t kernelHeaderLen = 18;
static const char kernelHeader[kernelHeaderLen+1] = "JDFTxCoulombKernel";

void CoulombKernelDesc::saveKernel(double* data, string filename) const
{
	FILE* fp = 0;
	if(filename.length())
	{	fp = fopen(filename.c_str(), "wb");
		if(!fp) logPrintf("WARNING: could not open %s for writing kernel.\n", filename.c_str());
		else logPrintf("Saving kernel to '%s' ... ", filename.c_str());
		logFlush();
	}
	
	//Write header:
	if(fp)
	{	fwrite(kernelHeader, sizeof(char), kernelHeaderLen, fp);
		fwrite(&R, sizeof(matrix3<>), 1, fp);
		fwrite(&S, sizeof(vector3<int>), 1, fp);
		fwrite(&isTruncated, sizeof(vector3<bool>), 1, fp);
		fwrite(&sigmaBorder, sizeof(vector3<>), 1, fp);
		fwrite(&omega, sizeof(double), 1, fp);
	}
	
	//Symmertize and write data:
	int nG = S[0]*S[1]*(1+S[2]/2);
	assert(nG > 0); //check for overflow
	std::vector<bool> done(nG, false); //mask for points whose symmetric image has already been output:
	std::vector<int> images(sym.size());
	
	size_t iStart=0, iStop=nG;
	THREAD_halfGspaceLoop
	(	if(!done[i])
		{	int nImages = 0;
			double dataSum = 0.;
			LOOP_symImagesG
			(	//process image point:
				dataSum += data[j];
				done[j] = true;
				images[nImages++] = j;
			)
			//Compute mean and symmetrize:
			double dataMean = dataSum / nImages;
			for(int iImage=0; iImage<nImages; iImage++)
				data[images[iImage]] = dataMean;
			//Output mean:
			if(fp) fwrite(&dataMean, sizeof(double), 1, fp);
		}
	)
	if(fp) { fclose(fp); logPrintf("Done.\n"); logFlush(); }
}

bool CoulombKernelDesc::loadKernel(double* data, string filename) const
{	if(!filename.size()) return false;
	FILE* fp = fopen(filename.c_str(), "rb");
	if(!fp)
	{	logPrintf("Could not open %s for reading kernel, it will be computed now.\n", filename.c_str());
		return false;
	}
	logPrintf("Reading kernel from %s ... ", filename.c_str()); logFlush();
	//Check header:
	char header[kernelHeaderLen];
	if(fread(header, sizeof(char), kernelHeaderLen, fp) != kernelHeaderLen)
	{	logPrintf("Failed: Premature EOF, kernel will be computed now.\n");
		fclose(fp); return false;
	}
	if(strncmp(header, kernelHeader, kernelHeaderLen))
	{	logPrintf("Failed: Invalid header, kernel will be computed now.\n");
		fclose(fp); return false;
	}
	//Check parameters:
	#define CHECKparam(paramType, param, paramName) \
		{	paramType paramIn; \
			if(fread(&paramIn, sizeof(param), 1, fp) != 1) \
			{	logPrintf("Failed: Premature EOF, kernel will be computed now.\n"); \
				fclose(fp); return false; \
			} \
			if(!(param == paramIn)) \
			{	logPrintf("Failed: mismatch in " paramName ", kernel will be computed now.\n"); \
				fclose(fp); return false; \
			} \
		}
	CHECKparam(matrix3<>, R, "lattice vectors")
	CHECKparam(vector3<int>, S, "sample count")
	CHECKparam(vector3<bool>, isTruncated, "truncated direction list")
	CHECKparam(vector3<>, sigmaBorder, "truncation border widths")
	CHECKparam(double, omega, "screening range parameter")
	#undef CHECKparam
	
	//Read data:
	size_t nG = S[0]*S[1]*(1+S[2]/2);
	std::vector<bool> done(nG, false); //mask for points whose symmetric image has already been output:
	
	size_t iStart=0, iStop=nG;
	THREAD_halfGspaceLoop
	(	if(!done[i])
		{	double dataMean;
			if(fread(&dataMean, sizeof(double), 1, fp) != 1)
			{	logPrintf("Failed: Premature EOF, kernel will be computed now.\n");
				fclose(fp); return false;
			}
			LOOP_symImagesG
			(	//save to image points:
				data[j] = dataMean;
				done[j] = true;
			)
		}
	)
	logPrintf("Done.\n"); logFlush();
	fclose(fp); return true;
}


void CoulombKernelDesc::computeIsolatedKernel(double* data, const WignerSeitz& ws) const
{	//Make sure all directions are truncated:
	for(int k=0; k<3; k++) assert(isTruncated[k]);
	//Count number of orthogonal pairs of directions:
	int nOrtho = 0;
	for(int k=0; k<3; k++)
		if(isOrthogonal(R.column(k), R.column((k+1)%2)))
			nOrtho++;
	//Chekc for speical case of cylinder (signaled by a negative borderSigma)
	for(int k=0; k<3; k++)
		if(sigmaBorder[k] < 0.)
		{	computeCylindrical(data, ws); return;
		}
	//Call appropriate specialized routine:
	switch(nOrtho)
	{	case 3: computeOrthorhombic(data, ws); return;
		case 2: computeMonoclinic(data, ws); return;
		default: computeTriclinic(data, ws);
	}
}


//--------- Triclinic geometry: no lattice vector orthogonal to other two ---------

namespace CoulombTriclinic
{
	//Set the fourier transform of the simplex theta function
	inline void simplexTheta_thread(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GT, double Vcell,
		complex* theta, const std::vector<Simplex<3>>* simplexArr, double sigma)
	{
		double volPrefac = 1./Vcell;
		THREAD_halfGspaceLoop
		(	vector3<> G = GT * iG; //reciprocal lattice vector in cartesian coords
			Simplex<3>::Point Gpoint({{ G[0], G[1], G[2] }}); //convert to Simplex<3>::Point
			double curTheta = 0.;
			for(const Simplex<3>& simplex: *simplexArr)
				curTheta += simplex.getTilde(Gpoint);
			theta[i] = volPrefac * curTheta * exp(-0.5*G.length_squared()*sigma*sigma);
		)
	}

	//Multiply by the erf-smoothened coulomb potential in real space, and initialize the index map
	inline void multErfCoulomb_thread(size_t iStart, size_t iStop, vector3<int> Spad, matrix3<> Rpad,
		double* data, const WignerSeitz* wsPad, double sigma, double dV,
		int* densePadMap, vector3<int> Sdense, const WignerSeitz* wsDense)
	{
		vector3<int> pitchPad;
		pitchPad[2] = 1;
		pitchPad[1] = pitchPad[2] * 2*(1+Spad[2]/2);
		pitchPad[0] = pitchPad[1] * Spad[1];
		vector3<int> pitchDense;
		pitchDense[2] = 1;
		pitchDense[1] = pitchDense[2] * 2*(1+Sdense[2]/2);
		pitchDense[0] = pitchDense[1] * Sdense[1];
		vector3<> invSpad, invSdense;
		for(int k=0; k<3; k++)
		{	invSpad[k] = 1./Spad[k];
			invSdense[k] = 1./Sdense[k];
		}
		matrix3<> hPad; //mesh offset vectors
		for(int k=0; k<3; k++)
			hPad.set_col(k, Rpad.column(k) / Spad[k]);
		matrix3<> hThPad = (~hPad) * hPad; //metric in mesh coordinates
		
		double a = sqrt(0.5)/sigma;
		const vector3<int>& S = Spad; //dimensions for THREAD_rLoop
		THREAD_rLoop
		(	//Compute index mappings:
			vector3<int> ivPad = wsPad->restrict(iv, Spad, invSpad); //position in mesh coordinates within padded WignerSeitz cell
			vector3<int> ivDense = wsDense->restrict(ivPad, Sdense, invSdense); //position in mesh coordinates within orig WS cell
			for(int k=0; k<3; k++)
			{	ivDense[k] = ivDense[k] % Sdense[k];
				if(ivDense[k]<0) ivDense[k] += Sdense[k];
			}
			int iPad = dot(iv, pitchPad);
			densePadMap[iPad] = dot(ivDense, pitchDense);
			//Multiply data by erf/r:
			double ar = a * sqrt(hThPad.metric_length_squared(ivPad)); //distance of minmal periodic image from origin
			data[iPad] *= (dV*a) * erf_by_x(ar); //include normalization factor for subsequent Fourier transform
		)
	}

	//Restrict from supercell to original cell, and add contribution from erfc/r
	inline void downSample_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
		const complex* in, const vector3<int>& Sdense, double* out, double sigma)
	{
		vector3<int> pitchDense;
		pitchDense[2] = 1;
		pitchDense[1] = pitchDense[2] * (1+Sdense[2]/2);
		pitchDense[0] = pitchDense[1] * Sdense[1];
		double hlfSigmaSq = 0.5*sigma*sigma;
		THREAD_halfGspaceLoop
		(	//Find corresponding point in dense cell fourier transform:
			int iDense = 0;
			for(int k=0; k<3; k++)
				iDense += pitchDense[k] * (iG[k]<0 ? iG[k]+Sdense[k] : iG[k]);
			//Store in smaller cell, along with short-ranged part:
			double Gsq = GGT.metric_length_squared(iG);
			out[i] = in[iDense].real() + (4*M_PI) * (Gsq ? (1.-exp(-hlfSigmaSq*Gsq))/Gsq : hlfSigmaSq);
		)
	}
}

void CoulombKernelDesc::computeTriclinic(double* data, const WignerSeitz& ws) const
{
	assert(sigmaBorder[0] == sigmaBorder[1]);
	assert(sigmaBorder[0] == sigmaBorder[2]);
	double borderWidth = sigmaBorder[0] * nSigmasPerWidth;
	double sigma = (ws.inRadius() - borderWidth) / nSigmasPerWidth;
	
	//Set up dense integration grids:
	logPrintf("Setting up FFT grids: ");
	double Gnyq = 10./sigmaBorder[0]; //lower bound on Nyquist frequency
	vector3<int> Sdense; //dense fft sample count
	matrix3<> Rpad; vector3<int> Spad; //padded lattice vectors and sample count
	for(int k=0; k<3; k++)
	{	Sdense[k] = std::max(S[k], 2*int(ceil(Gnyq * R.column(k).length() / (2*M_PI))));
		while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
		//Pad the super cell by border width:
		Spad[k] = Sdense[k] + 2*ceil(Sdense[k]*borderWidth/R.column(k).length());
		while(!fftSuitable(Spad[k])) Spad[k]+=2;
		Rpad.set_col(k, R.column(k) * (Spad[k]*1./Sdense[k]));
	}
	logPrintf("[ %d %d %d ], and [ %d %d %d ] padded.\n",
		Sdense[0], Sdense[1], Sdense[2], Spad[0], Spad[1], Spad[2]);
	size_t nG = S[0] * (S[1] * size_t(1 + S[2]/2));
	size_t nGdense = Sdense[0] * (Sdense[1] * size_t(1+Sdense[2]/2));
	size_t nGpad = Spad[0] * (Spad[1] * size_t(1+Spad[2]/2)); //number of symmetry reduced reciprocal lattice vectors
	size_t nrPad = Spad[0] * (Spad[1] * size_t(Spad[2])); //number of real space points
	double detRpad = fabs(det(Rpad));
	
	//Construct padded Wigner-Seitz cell, and check border:
	logPrintf("For padded lattice, "); WignerSeitz wsPad(Rpad);
	std::vector<vector3<>> vArr = ws.getVertices();
	matrix3<> invRpad = inv(Rpad);
	for(vector3<> v: vArr)
		if(wsPad.boundaryDistance(wsPad.restrict(invRpad*v)) < 0.9*borderWidth)
			die("\nPadded Wigner-Seitz cell does not fit inside Wigner-Seitz cell of padded lattice.\n"
				"This can happen for reducible lattice vectors; HINT: the reduced lattice vectors,\n"
				"if different from the input lattice vectors, are printed during Symmetry setup.\n");
	
	//Plan Fourier transforms:
	complex* padArr = (complex*)fftw_malloc(sizeof(complex)*nGpad);
	complex* denseArr = (complex*)fftw_malloc(sizeof(complex)*nGdense);
	int* densePadMap = new int[nGpad*2]; //index map from padded to dense grid
	double* padRealArr = (double*)padArr;
	double* denseRealArr = (double*)denseArr;
	if(!padArr || !denseArr || !densePadMap)
		die("Insufficient memory (need %.1fGB). Hint: try increasing border width.\n",
			(nGpad+nGdense)*1e-9*sizeof(complex) + nGpad*1e-9*2*sizeof(int));
	logPrintf("Planning fourier transforms ... "); logFlush();
	fftw_plan_with_nthreads(nProcsAvailable);
	fftw_plan fftPlanC2R = fftw_plan_dft_c2r_3d(Spad[0], Spad[1], Spad[2], (fftw_complex*)padArr, padRealArr, FFTW_ESTIMATE);
	fftw_plan fftPlanR2C = fftw_plan_dft_r2c_3d(Sdense[0], Sdense[1], Sdense[2], denseRealArr, (fftw_complex*)denseArr, FFTW_ESTIMATE);
	logPrintf("Done.\n");
	
	//Initialize smoothed theta function on padded lattice:
	logPrintf("Computing truncation shape function ... "); logFlush();
	std::vector<Simplex<3>> sArr = ws.getSimplices();
	threadLaunch(CoulombTriclinic::simplexTheta_thread, nGpad, Spad, (2.*M_PI)*(~invRpad), detRpad, padArr, &sArr, sigmaBorder[0]);
	fftw_execute(fftPlanC2R);
	logPrintf("Done.\n");
	
	//Multiply by long-ranged erf/r and initialize index map (periodic map via Wigner-Seitz restrictions):
	logPrintf("Applying truncation to Coulomb kernel ... "); logFlush();
	threadLaunch(CoulombTriclinic::multErfCoulomb_thread, nrPad, Spad, Rpad, padRealArr, &wsPad, sigma, detRpad/nrPad, densePadMap, Sdense, &ws);
	//Collect values from points in padded grid into dense grid using index map:
	//(This involves scatter operations, so cannot be threaded (easily))
	memset(denseArr, 0, sizeof(complex)*nGdense);
	int pitchPad2 = 2*(1+Spad[2]/2);
	for(int i01=0; i01<Spad[0]*Spad[1]; i01++)
	{	int iPad = i01*pitchPad2;
		for(int i2=0; i2<Spad[2]; i2++)
		{	denseRealArr[densePadMap[iPad]] += padRealArr[iPad];
			iPad++;
		}
	}
	fftw_free(padArr);
	delete[] densePadMap;
	fftw_execute(fftPlanR2C);
	
	//Restrict to original sample count and add short-ranged erfc/r in reciprocal space:
	matrix3<> G = (2.*M_PI) * inv(R);
	matrix3<> GGT = G * (~G);
	threadLaunch(CoulombTriclinic::downSample_thread, nG, S, GGT, denseArr, Sdense, data, sigma);
	fftw_free(denseArr);
	fftw_destroy_plan(fftPlanC2R);
	fftw_destroy_plan(fftPlanR2C);
	logPrintf("Done.\n");
}



//--------- Monoclinic geometry: one lattice vector orthogonal to other two ---------
void CoulombKernelDesc::computeMonoclinic(double* data, const WignerSeitz& ws) const
{	//TODO: Optimize special case
	computeTriclinic(data, ws);
}



//--------- Orthorhombic geometry: all three lattice vectors orthogonal ---------
void CoulombKernelDesc::computeOrthorhombic(double* data, const WignerSeitz& ws) const
{	//TODO: Optimize special case
	computeTriclinic(data, ws);
}



//--------- Cylindrical geometry: truncation on a finite length cylinder  ---------
void CoulombKernelDesc::computeCylindrical(double* data, const WignerSeitz& ws) const
{	
	//Determine and verify geometry:
	int iDir = -1;
	for(int k=0; k<3; k++)
		if(sigmaBorder[k] > 0.)
			iDir = k;
	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	double Rc = -sigmaBorder[jDir];
	assert(Rc > 0.);
	assert(Rc <= ws.inRadius(iDir));
	assert(Rc == -sigmaBorder[kDir]);
	assert(isOrthogonal(R.column(iDir), R.column(jDir)));
	assert(isOrthogonal(R.column(iDir), R.column(kDir)));
	
	die("Not yet implemented.\n");
}

