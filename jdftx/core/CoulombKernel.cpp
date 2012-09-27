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

//Number of gaussian sigma's to include in border width,
//so as to make overlap negligible at working precision
const double CoulombKernelDesc::nSigmasPerWidth = 10.;

double CoulombKernelDesc::getMaxSigma(double L, double sigmaOther)
{	double sigmaHypotMax = L / nSigmasPerWidth;
	double maxSigmaSq = sigmaHypotMax*sigmaHypotMax - sigmaOther*sigmaOther;
	return maxSigmaSq>0 ? sqrt(maxSigmaSq) : 0.;
}

//getMaxSigma with an overlap check for the exchange kernel cases:
double getMaxSigma_overlapCheck(double L, double sigmaOther)
{	double sigma = CoulombKernelDesc::getMaxSigma(L, sigmaOther);
	if(sigma < sigmaOther)
	{	sigma = sigmaOther;
		const double overlapWarnThreshold = 1e-13;
		const double sigmaMax = 0.5*L/sqrt(log(1./overlapWarnThreshold));
		if(sigma >= sigmaMax)
		{	double overlapEstimate = exp(-0.25*pow(L/sigma,2));
			logPrintf(
				"WARNING: Gaussian border and short-ranged part of potential overlap ~ %.0le\n"
				"         This is most likely important only to the paranoid; this warning\n"
				"         will be suppressed for sigma < %lg for the current geometry.\n",
				overlapEstimate, sigmaMax);
		}
	}
	return sigma;
}


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
	std::vector<matrix3<int>> symLattice = ::getSymmetries(R); //symmetries of the Bravais lattice
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

std::vector<matrix3<int>> CoulombKernelDesc::getSymmetries() const
{	return sym;
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


void CoulombKernelDesc::computeKernel(double* data, const WignerSeitz& ws, string filename) const
{	//Count number of orthogonal pairs of directions:
	int nOrtho = 0;
	for(int k=0; k<3; k++)
		if(isOrthogonal(R.column(k), R.column((k+1)%2)))
			nOrtho++;
	//Try loading from file:
	if(loadKernel(data, filename))
		return; //successfully loaded
	//Call appropriate specialized routine:
	if(nOrtho<2) { computeNonOrtho(data, ws); return; }
	else { computeRightPrism(data, ws); return; }
	//Save kernel to file:
	saveKernel(data, filename);
}

//! Compute erfc(omega r)/r - erfc(a r)/r
inline double screenedPotential(double r, double a, double omega)
{	double result = a * erf_by_x(a * r);
	if(omega) result -= omega * erf_by_x(omega * r);
	return result;
}

//--------- General case: no lattice vector orthogonal to other two ---------

namespace CoulombNonOrtho
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
		double* data, const WignerSeitz* wsPad, double sigma, double omega, double dV,
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
			//Multiply data by long-range part of optionally screened potential:
			double r = sqrt(hThPad.metric_length_squared(ivPad)); //distance of minmal periodic image from origin
			data[iPad] *= dV * screenedPotential(r, a, omega); //include normalization for subsequent Fourier transform
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

void CoulombKernelDesc::computeNonOrtho(double* data, const WignerSeitz& ws) const
{
	for(int k=0; k<3; k++) assert(isTruncated[k]); //Make sure all directions are truncated
	assert(sigmaBorder[0] == sigmaBorder[1]);
	assert(sigmaBorder[0] == sigmaBorder[2]);
	double borderWidth = sigmaBorder[0] * nSigmasPerWidth;
	double sigma = getMaxSigma_overlapCheck(ws.inRadius(), sigmaBorder[0]);
	logPrintf("Using arbitrary geometry algorithm (3D simplices).\n");
	
	//Set up dense integration grids:
	logPrintf("Setting up FFT grids: ");
	double Gnyq = nSigmasPerWidth/sigmaBorder[0]; //lower bound on Nyquist frequency
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
	threadLaunch(CoulombNonOrtho::simplexTheta_thread, nGpad, Spad, (2.*M_PI)*(~invRpad), detRpad, padArr, &sArr, sigmaBorder[0]);
	fftw_execute(fftPlanC2R);
	logPrintf("Done.\n");
	
	//Multiply by long-ranged erf/r and initialize index map (periodic map via Wigner-Seitz restrictions):
	logPrintf("Applying truncation to Coulomb kernel ... "); logFlush();
	threadLaunch(CoulombNonOrtho::multErfCoulomb_thread, nrPad, Spad, Rpad, padRealArr, &wsPad, sigma, omega, detRpad/nrPad, densePadMap, Sdense, &ws);
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
	threadLaunch(CoulombNonOrtho::downSample_thread, nG, S, GGT, denseArr, Sdense, data, sigma);
	fftw_free(denseArr);
	fftw_destroy_plan(fftPlanC2R);
	fftw_destroy_plan(fftPlanR2C);
	logPrintf("Done.\n");
}

//--------- Right-prism geometry: at least one lattice vector orthogonal to other two ---------

//1D fourier transform of sharply truncated erf/r:
struct TruncatedErfTilde
{
	TruncatedErfTilde(double k, double sigma, double omega, double hlfL)
	: aPot(sqrt(0.5)/sigma), omega(omega), absTol(1e-13), relTol(1e-12*std::max(0.1,k*hlfL))
	{	iWS = gsl_integration_workspace_alloc(nBisections);
		qawoTable = gsl_integration_qawo_table_alloc(k, hlfL, GSL_INTEG_COSINE, nBisections);
	}
	~TruncatedErfTilde()
	{	gsl_integration_qawo_table_free(qawoTable);
		gsl_integration_workspace_free(iWS);
	}
	double operator()(double rho)
	{	double params[3] = { rho, aPot, omega }, result, err;
		gsl_function f = { &integrand, params };
		gsl_integration_qawo(&f, 0., absTol, relTol, nBisections, iWS, qawoTable, &result, &err);
		return 2. * result; //factor of 2 for symmetric z<0
	}
private:
	const double aPot, omega, absTol, relTol;
	static const int nBisections = 50;
	gsl_integration_workspace* iWS; //integration workspace
	gsl_integration_qawo_table* qawoTable; //integration table for oscillatory integrand
	
	static double integrand(double z, void* params)
	{	const double* p = (const double*)params;
		const double& rho     = p[0]; //distance from axis
		const double& aPot    = p[1]; //sqrt(0.5)/sigma
		const double& omega   = p[2]; //erf-screening parameter in potential (0 for unscreened)
		double r = hypot(z, rho);
		return screenedPotential(r, aPot, omega);
	}
};

//Threaded computation of kernel for RightPrism geometry
struct CoulombRightPrism
{	
	//Data arrays:
	complex* padArr; double* padRealArr; //fft array on padded grid
	complex* denseArr; double* denseRealArr; //fft array on dense (unpadded) grid
	double* Vc; //output kernel (3D fftw c2r layout)
	fftw_plan fftPlanC2R, fftPlanR2C;
	const double* radialPot; double drho; int nSamples; //Radial potential look-up table
	
	//Geometry:
	int iDir, jDir, kDir;
	vector3<int> S, Spad, Sdense;
	matrix3<> Rplanar, RpadPlanar, GGT;
	std::vector<Simplex<2>>* simplexArr; //2D Wigner-Seitz cell simplicial tesselation
	const WignerSeitz *wsDense, *wsPad; //Wigner-Seitz cells
	double sigma, omega; vector3<> sigmaBorder; vector3<bool> isTruncated;
	
	
	void computePlane(int iPlane)
	{
		//Initialize look-up table (in rho) for the axial Fourier transform of truncated potential:
		std::vector<double> coeff; double drhoInv = 1./drho; bool logSpline = false;
		{	std::vector<double> samples(nSamples);
			if(radialPot) //smooth truncation along iDir
				samples.assign(radialPot+iPlane*nSamples, radialPot+(iPlane+1)*nSamples);
			else //sharp truncation or no truncation
			{	double hlfL = 0.5*Rplanar.column(iDir).length();
				double kz = iPlane * (M_PI/hlfL);
				if(isTruncated[iDir]) //sharp truncation
				{	TruncatedErfTilde erfTilde(kz, sigma, omega, hlfL);
					for(int i=0; i<nSamples; i++)
						samples[i] = erfTilde(i*drho);
				}
				else //no truncation
				{	logSpline = (kz != 0.); //Logarithmic spline more stable for large kz
					Cbar cbar;
					for(int i=0; i<nSamples; i++)
					{	double s = cbar(kz, sigma, i*drho);
						samples[i] = logSpline ? log(s) : s;
					}
				}
			}
			//Initialize coefficients:
			coeff = QuinticSpline::getCoeff(samples);
		}
		
		//Compute smoothed WignerSeitz-shaped theta function in Fourier space (non-orthorhombic case):
		if(padArr) //exclude orthorhombic case
		{	const double invApad = RpadPlanar.column(iDir).length() / fabs(det(RpadPlanar));
			const double sigmaBorderSq = pow(sigmaBorder[jDir], 2);
			assert(sigmaBorder[jDir] == sigmaBorder[kDir]);
			vector3<int> iGpad(0,0,0);
			matrix3<> GTpadPlanar = (2*M_PI) * ~inv(RpadPlanar);
			complex* theta = padArr;
			for(iGpad[jDir]=0;;)
			{	for(iGpad[kDir]=0; iGpad[kDir]<=Spad[kDir]/2; iGpad[kDir]++)
				{	vector3<> G = GTpadPlanar * iGpad; //reciprocal lattice vector in planar cartesian coords
					Simplex<2>::Point Gpoint({{ G[0], G[1] }}); //convert to Simplex<2>::Point (third direction is truncated)
					double curTheta = 0.;
					for(const Simplex<2>& simplex: *simplexArr)
						curTheta += simplex.getTilde(Gpoint);
					*(theta++) = invApad * curTheta * exp(-0.5*G.length_squared()*sigmaBorderSq);
				}
				iGpad[jDir]++;
				if(2*iGpad[jDir]>Spad[jDir]) iGpad[jDir] -= Spad[jDir];
				if(iGpad[jDir]==0) break;
			}
			fftw_execute_dft_c2r(fftPlanC2R, (fftw_complex*)padArr, padRealArr);
		}
		
		//Multiply by 1D fourier transform of potential in real-space and fold into dense array:
		if(denseArr) //exclude capped-cylinder case
		{	int pitchPad = 2*(1+Spad[kDir]/2);
			int pitchDense = 2*(1+Sdense[kDir]/2);
			vector3<> invSpad, invSdense;
			for(int k=0; k<3; k++)
			{	invSpad[k] = 1./Spad[k];
				invSdense[k] = 1./Sdense[k];
			}
			matrix3<> hPad; //mesh offset vectors
			vector3<> haBorder; //sqrt(0.5)/sigmaBorder in mesh coordinates
			for(int k=0; k<3; k++)
			{	hPad.set_col(k, RpadPlanar.column(k) / Spad[k]);
				haBorder[k] = hPad.column(k).length() * sqrt(0.5)/sigmaBorder[k];
			}
			#define SmoothTheta1D(k) 0.5*erfc(haBorder[k] * (abs(ivPad[k])-Sdense[k]/2))
			matrix3<> hThPad = (~hPad) * hPad; //metric in mesh coordinates
			double dA = fabs(det(Rplanar)) / (Rplanar.column(iDir).length() * Sdense[jDir] * Sdense[kDir]);
			vector3<int> iv(0,0,0);
			memset(denseArr, 0, sizeof(complex)*Sdense[jDir]*(1+Sdense[kDir]/2));
			for(iv[jDir]=0; iv[jDir]<Spad[jDir]; iv[jDir]++)
				for(iv[kDir]=0; iv[kDir]<Spad[kDir]; iv[kDir]++)
				{	//Compute index mappings:
					vector3<int> ivPad = wsPad->restrict(iv, Spad, invSpad); //position in mesh coordinates within padded WignerSeitz cell
					vector3<int> ivDense = wsDense->restrict(ivPad, Sdense, invSdense); //position in mesh coordinates within orig WS cell
					for(int k=0; k<3; k++)
					{	ivDense[k] = ivDense[k] % Sdense[k];
						if(ivDense[k]<0) ivDense[k] += Sdense[k];
					}
					int iPad = iv[kDir] + pitchPad * iv[jDir];
					int iDense = ivDense[kDir] + pitchDense * ivDense[jDir];
					//Accumulate theta multiplied by potential to the mapped position:
					double rho = sqrt(hThPad.metric_length_squared(ivPad)); //distance of minmal periodic image from origin
					double splineVal = QuinticSpline::value(coeff.data(), drhoInv*rho); //evaluate radial spline
					double pot = logSpline ? exp(splineVal) : splineVal; //1D fourier transform of 1D-truncated potential
					double theta = padArr //2D theta function
						? padRealArr[iPad] //From simplices in fourier space
						: SmoothTheta1D(jDir) * SmoothTheta1D(kDir); //Analytical (orthorhombic case)
					denseRealArr[iDense] += dA * pot * theta;
				}
			#undef SmoothTheta1D
			fftw_execute_dft_r2c(fftPlanR2C, denseRealArr, (fftw_complex*)denseArr);
		}
		
		//Add analytic short-ranged parts in fourier space (and down-sample to final resolution):
		int pitchDense = 1+Sdense[kDir]/2;
		vector3<int> pitch;
		pitch[2] = 1;
		pitch[1] = pitch[2] * (1 + S[2]/2);
		pitch[0] = pitch[1] * S[1];
		vector3<int> iG; iG[iDir] = iPlane;
		for(iG[jDir]=1-S[jDir]/2; iG[jDir]<=S[jDir]/2; iG[jDir]++)
			for(iG[kDir]=0; iG[kDir]<=S[kDir]/2; iG[kDir]++)
			{	//Collect the data from the dense grid transform:
				int iDense = iG[kDir] + pitchDense*(iG[jDir]<0 ? iG[jDir]+Sdense[jDir] : iG[jDir]);
				double curV = denseArr[iDense].real();
				//Add the analytical short-ranged part:
				double Gsq = GGT.metric_length_squared(iG);
				curV += (4*M_PI) * (Gsq ? (1.-exp(-0.5*sigma*sigma*Gsq))/Gsq : 0.5*sigma*sigma);
				//Save to the appropriate locations in Vc (4 sign combinations):
				for(int si=0; si<2; si++)
				{	for(int sk=0; sk<2; sk++)
					{	vector3<int> iv = iG;
						if(si && iv[iDir]) { iv[iDir] = S[iDir] - iv[iDir]; }
						if(sk && iv[kDir]) { iv[kDir] = S[kDir] - iv[kDir]; iv[jDir] = -iv[jDir]; }
						if(iv[jDir] < 0) iv[jDir] += S[jDir];
						if(iv[2] <= S[2]/2)
							Vc[dot(iv, pitch)] = curV;
					}
				}
			}
	}
	
	static void thread(int iThread, int nThreads, CoulombRightPrism* crpArr,
		int nPlanes, int* nPlanesDone, std::mutex* m)
	{
		while(true)
		{	//Get next available job:
			m->lock();
			int iPlane = (*nPlanesDone)++;
			m->unlock();
			if(iPlane >= nPlanes)
				break; //job queue empty
			//Perform job:
			crpArr[iThread].computePlane(iPlane);
		}
	}
};

void CoulombKernelDesc::computeRightPrism(double* data, const WignerSeitz& ws) const
{	//This function has 3 distinct modes which share common pieces:
	// * 2D-truncated wire geometry
	// * Capped cylinder geometry
	// * 3D truncated geometry with at least one axis orthogonal
	//   to the other two, which in turn has two sub-modes:
	//    - Orthorhombic (all three axes orthogonal)
	//    - Hexagonal right-prism (general 2D wigner-sitz cell is a hexagon)

	//Determine mode and pick axis:
	int iDir = -1; bool cylinderMode = false;
	double sigmaBorderMin = *std::min_element(&sigmaBorder[0], &sigmaBorder[0]+3);
	double sigmaBorderMax = *std::max_element(&sigmaBorder[0], &sigmaBorder[0]+3);
	//--- Check for wire geometry:
	for(int k=0; k<3; k++)
		if(!isTruncated[k])
		{	assert(iDir < 0); //there should be at most one untruncated direction
			assert(!omega); //not implemented for unregularized screened-exchange
			iDir = k; //axis is along untruncated direction
		}
	//--- Check for capped cylinder:
	if(iDir<0 && sigmaBorderMin<0.)
	{	for(int k=0; k<3; k++)
			if(sigmaBorder[k] >= 0.)
				iDir = k; //axis is along direction with non-negative sigma
		assert(iDir >= 0);
		cylinderMode = true;
		die("Cylinder geometry exchange-regularization implementation incomplete.\n");
	}
	//--- Check for finite prism:
	if(iDir<0)
	{	//Pick iDir to be direction orthogonal to ther two with smallest sigmaBorder:
		double sigmaAxis = DBL_MAX;
		for(int k=0; k<3; k++)
			if( isOrthogonal(R.column(k), R.column((k+1)%3))
			&&  isOrthogonal(R.column(k), R.column((k+2)%3))
			&&  (sigmaBorder[k] < sigmaAxis) )
			{	iDir = k;
				sigmaAxis = sigmaBorder[k];
			}
		assert(iDir >= 0);
		assert(sigmaAxis >= 0);
	}
	
	//Check transverse geometry:
	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	assert(isOrthogonal(R.column(iDir), R.column(jDir)));
	assert(isOrthogonal(R.column(iDir), R.column(kDir)));
	double sigmaPlaneMin = std::min(sigmaBorder[jDir], sigmaBorder[kDir]);
	assert(sigmaPlaneMin > 0.); //this determines Nyquist frequency for planar FFT
	bool jkOrtho = isOrthogonal(R.column(jDir), R.column(kDir));
	if(!jkOrtho || cylinderMode) assert(sigmaBorder[jDir] == sigmaBorder[kDir]);
	if(cylinderMode) assert(-sigmaBorder[jDir] <= ws.inRadius(iDir));
	double sigma = getMaxSigma_overlapCheck(ws.inRadius(), sigmaBorderMax);
	
	//Print mode-dependent startup message:
	string dirName(3,'0'); dirName[iDir]='1';
	if(cylinderMode)
		logPrintf("Computing kernel truncated on finite cylinder (%s-axis x %lg bohrs radius)\n",
			dirName.c_str(), -sigmaBorder[jDir]);
	else
	{	string axisMode = isTruncated[iDir] ? "truncated" : "infinite";
		string geomName = jkOrtho ? "Orthorhombic" : "Right-Prism";
		string geomDesc = jkOrtho ? "erfc x erfc" : "2D simplices";
		logPrintf("Using %s geometry algorithm (%s %s-axis x %s)\n",
			geomName.c_str(), axisMode.c_str(), dirName.c_str(), geomDesc.c_str());
	}
	
	//Set up dense integration grids:
	vector3<int> Sdense; //dense fft sample count
	matrix3<> Rpad; vector3<int> Spad; //padded lattice vectors and sample count
	if(cylinderMode) 
	{	//Don't need any 2D/3D grids besides the original one:
		Sdense = S;
		Spad = S;
		Rpad = R;
	}
	else
	{	logPrintf("Setting up FFT grids: ");
		double Gnyq = nSigmasPerWidth/sigmaPlaneMin; //lower bound on Nyquist frequency
		for(int k=0; k<3; k++)
			if(k == iDir)
			{	Sdense[k] = S[k];
				Spad[k] = S[k];
				Rpad.set_col(k, R.column(k));
			}
			else
			{	Sdense[k] = std::max(S[k], 2*int(ceil(Gnyq * R.column(k).length() / (2*M_PI))));
				while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
				//Pad the super cell by border width:
				double borderWidth = sigmaBorder[k] * nSigmasPerWidth;
				Spad[k] = Sdense[k] + 2*ceil(Sdense[k]*borderWidth/R.column(k).length());
				while(!fftSuitable(Spad[k])) Spad[k]+=2;
				Rpad.set_col(k, R.column(k) * (Spad[k]*1./Sdense[k]));
			}
		logPrintf("%d x %d, and %d x %d padded.\n", Sdense[jDir], Sdense[kDir], Spad[jDir], Spad[kDir]);
	}
	int nGdense = Sdense[jDir] * (1+Sdense[kDir]/2);
	int nGpad = Spad[jDir] * (1+Spad[kDir]/2); //number of symmetry reduced planar reciprocal lattice vectors
	int nrPad = Spad[jDir] * Spad[kDir]; //number of planar real space points
	matrix3<> G = (2*M_PI)*inv(R);
	matrix3<> GGT = G * (~G);
	
	//Construct padded Wigner-Seitz cell, and check border:
	WignerSeitz* wsPad = 0;
	if(!cylinderMode)
	{	logPrintf("For padded lattice, "); wsPad = new WignerSeitz(Rpad);
		double borderWidthMin = sigmaPlaneMin * nSigmasPerWidth;
		std::vector<vector3<>> vArr = ws.getVertices();
		matrix3<> invRpad = inv(Rpad);
		for(vector3<> v: vArr)
			if(wsPad->boundaryDistance(wsPad->restrict(invRpad*v), iDir) < 0.9*borderWidthMin)
				die("\nPadded Wigner-Seitz cell does not fit inside Wigner-Seitz cell of padded lattice.\n"
					"This can happen for reducible lattice vectors; HINT: the reduced lattice vectors,\n"
					"if different from the input lattice vectors, are printed during Symmetry setup.\n");
	}
	//Plan Fourier transforms:
	fftw_plan fftPlanC2R = 0, fftPlanR2C = 0;
	if(!cylinderMode)
	{	assert(nrPad > 0); //overflow check
		complex* tempArr = (complex*)fftw_malloc(sizeof(complex)*nGpad);
		#define INSUFFICIENT_MEMORY_ERROR \
			die("Insufficient memory (will need %.1fGB). Hint: try increasing border width.\n", \
				((jkOrtho ? 0 : nGpad) + nGdense)*1e-9*nProcsAvailable*sizeof(complex));
		if(!tempArr) INSUFFICIENT_MEMORY_ERROR
		logPrintf("Planning fourier transforms ... "); logFlush();
		fftw_plan_with_nthreads(1); //Multiple simultaneous single threaded fourier transforms
		fftPlanC2R = fftw_plan_dft_c2r_2d(Spad[jDir], Spad[kDir], (fftw_complex*)tempArr, (double*)tempArr, FFTW_ESTIMATE);
		fftPlanR2C = fftw_plan_dft_r2c_2d(Sdense[jDir], Sdense[kDir], (double*)tempArr, (fftw_complex*)tempArr, FFTW_ESTIMATE);
		fftw_free(tempArr);
		logPrintf("Done.\n");
	}

	//Compute radial look-up table from a 1D FFT if iDir is truncated smoothly:
	double drho = 0.03 * sigma;
	double rhoMax = cylinderMode ? -sigmaBorder[jDir] : wsPad->circumRadius(iDir);
	int nSamples = ceil(rhoMax/drho) + 10;
	double* radialPot = 0;
	if(isTruncated[iDir] && sigmaBorder[iDir] > 0.)
	{	radialPot = new double[nSamples * (1+S[iDir]/2)];
		if(!radialPot) INSUFFICIENT_MEMORY_ERROR
		//Setup a fine grid for z-integration:
		double Gnyq = nSigmasPerWidth/std::min(sigma,sigmaBorder[iDir]);
		double hlfL = 0.5 * R.column(iDir).length(); //half-length along axis
		int nAxis = std::max(S[iDir]/2, int(ceil(Gnyq * hlfL / M_PI)));
		while(!fftSuitable(2*nAxis)) nAxis++;
		double hAxis = hlfL / nAxis;
		int nPad = std::min(int(sigmaBorder[iDir]*nSigmasPerWidth/hAxis), nAxis);
		//Plan FFT (reduces to Discrete Cosine Transform):
		double* fftArr = (double*)fftw_malloc(sizeof(double)*(nAxis+1));
		if(!fftArr) INSUFFICIENT_MEMORY_ERROR
		fftw_plan fftPlanDCT1 = fftw_plan_r2r_1d(nAxis+1, fftArr, fftArr, FFTW_REDFT00, FFTW_ESTIMATE);
		//1D FFT for each radial location:
		double aPot = sqrt(0.5)/sigma;
		double aBorder = sqrt(0.5)/sigmaBorder[iDir];
		for(int iRho=0; iRho<nSamples; iRho++)
		{	double rho = drho * iRho;
			//Fill fftArr in real space:
			for(int iz=0; iz<nAxis+nPad; iz++)
			{	double z = iz * hAxis;
				double r = hypot(z, rho);
				double term = hAxis //sum -> integral
					* screenedPotential(r, aPot, omega) //long-ranged part of potential (smooth)
					* (0.5 * erfc(aBorder*(z-hlfL))); //smooth theta function for truncating iDir
				if(iz<=nAxis) fftArr[iz] = term;
				else fftArr[2*nAxis-iz] += term; //fold tail using periodicity
			}
			fftArr[nAxis] *= 2; //this sample gets half-weighted by the FFT routine
			fftw_execute(fftPlanDCT1);
			//Save 1D fourier transform to radialPot:
			for(int iz=0; iz<=S[iDir]/2; iz++)
				radialPot[iRho+iz*nSamples] = fftArr[iz];
		}
		//Cleanup:
		fftw_destroy_plan(fftPlanDCT1);
		fftw_free(fftArr);
	}
	
	//Setup threads for initializing each plane perpendicular to axis:
	logPrintf("Computing truncated kernel ... "); logFlush();
	std::vector<CoulombRightPrism> crpArr(nProcsAvailable);
	std::vector<Simplex<2>> simplexArr = ws.getSimplices(iDir);
	for(CoulombRightPrism& c: crpArr)
	{	//Allocate data arrays (padded array not needed for orthorhombic case):
		if(cylinderMode)
		{	c.padArr = c.denseArr = 0;
		}
		else
		{	c.padArr = jkOrtho ? 0 : (complex*)fftw_malloc(sizeof(complex)*nGpad);
			c.denseArr = (complex*)fftw_malloc(sizeof(complex)*nGdense);
			if(!c.denseArr) INSUFFICIENT_MEMORY_ERROR
			#undef INSUFFICIENT_MEMORY_ERROR
		}
		c.padRealArr = (double*)c.padArr;
		c.denseRealArr = (double*)c.denseArr;
		c.Vc = data;
		c.fftPlanC2R = fftPlanC2R;
		c.fftPlanR2C = fftPlanR2C;
		c.radialPot = radialPot; c.drho = drho; c.nSamples = nSamples;
		//Copy geometry definitions:
		c.iDir = iDir; c.jDir = jDir; c.kDir = kDir;
		c.S = S; c.Sdense = Sdense; c.Spad = Spad;
		c.Rplanar = ws.getRplanar(iDir);
		c.RpadPlanar = (wsPad ? wsPad->getRplanar(iDir) : c.Rplanar);
		c.GGT = GGT;
		c.simplexArr = &simplexArr;
		c.wsDense = &ws; c.wsPad = wsPad;
		c.sigma = sigma; c.omega = omega;
		c.sigmaBorder = sigmaBorder; c.isTruncated = isTruncated;
	}
	
	//Launch threads:
	std::mutex mJobCount; int nPlanesDone = 0; //for job management
	threadLaunch(CoulombRightPrism::thread, 0, crpArr.data(), 1+S[iDir]/2, &nPlanesDone, &mJobCount);
	if(fftPlanC2R) fftw_destroy_plan(fftPlanC2R);
	if(fftPlanR2C) fftw_destroy_plan(fftPlanR2C);
	if(radialPot) delete[] radialPot;
	if(wsPad) delete wsPad;
	
	//Cleanup threads:
	for(CoulombRightPrism& c: crpArr)
	{	if(!cylinderMode)
		{	if(!jkOrtho) fftw_free(c.padArr);
			fftw_free(c.denseArr);
		}
	}
	logPrintf("Done.\n");
}
