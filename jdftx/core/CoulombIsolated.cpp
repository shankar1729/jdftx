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

#include <core/CoulombIsolated.h>
#include <core/Coulomb_internal.h>
#include <core/Operators.h>
#include <core/Util.h>
#include <core/Thread.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>

//! Analog of ewald sum for isolated systems
//! (no Ewald trick required, just for consistent naming)
struct EwaldIsolated
{	const GridInfo& gInfo;
	const WignerSeitz& ws;
	bool wsTruncated; //true => Wigner-Seitz truncation, false => spherical
	double criticalDist; //borderWidth for Wigner-Seitz, Rc for spherical
	
	EwaldIsolated(const GridInfo& gInfo, const WignerSeitz& ws, bool wsTruncated, double criticalDist)
	: gInfo(gInfo), ws(ws), wsTruncated(wsTruncated), criticalDist(criticalDist)
	{
	}
	
	double energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
	{	if(!pointCharges.size()) return 0.;
		double E = 0.;
		//Shift all points into a Wigner-Seitz cell centered on one of the atoms; choice of this atom
		//is irrelevant if every atom lies in the WS cell of the other with a consistent translation:
		vector3<> pos0 = pointCharges[0].pos;
		for(Coulomb::PointCharge& pc: pointCharges)
			pc.pos = pos0 + ws.restrict(pc.pos - pos0);
		//Loop over all pairs of pointcharges:
		for(unsigned i=0; i<pointCharges.size(); i++)
		{	Coulomb::PointCharge& pc1 = pointCharges[i];
			for(unsigned j=0; j<i; j++)
			{	Coulomb::PointCharge& pc2 = pointCharges[j];
				vector3<> x = pc1.pos - pc2.pos; //lattice coords
				double rSq = gInfo.RTR.metric_length_squared(x), r = sqrt(rSq);
				if(wsTruncated)
				{	if(ws.boundaryDistance(x) <= criticalDist)
						die("Separation between atoms %d and %d lies in the truncation border + margin.\n", i, j);
				}
				else
				{	if(r >= criticalDist)
						die("Atoms %d and %d are separated by r = %lg >= Rc-ionMargin = %lg bohrs.\n", i, j, r, criticalDist);
				}
				double dE = (pc1.Z * pc2.Z) / r;
				vector3<> dF = (gInfo.RTR * x) * (dE/rSq);
				E += dE;
				pc1.force += dF;
				pc2.force -= dF;
			}
		}
		return E;
	}
};

//----------------- class CoulombIsolated ---------------------

//Set the fourier transform of the simplex theta function
inline void simplexTheta_thread(int iStart, int iStop, const vector3<int> S, const matrix3<> GT, double Vcell,
	fftw_complex* theta, const std::vector<Simplex<3>>* simplexArr, double sigma)
{
	double volPrefac = 1./Vcell;
	THREAD_halfGspaceLoop
	(	vector3<> G = GT * iG; //reciprocal lattice vector in cartesian coords
		Simplex<3>::Point Gpoint({{ G[0], G[1], G[2] }}); //convert to Simplex<3>::Point
		double curTheta = 0.;
		for(const Simplex<3>& simplex: *simplexArr)
			curTheta += simplex.getTilde(Gpoint);
		theta[i][0] = volPrefac * curTheta * exp(-0.5*G.length_squared()*sigma*sigma);
		theta[i][1] = 0.; //no imaginary part by inversion symmetry
	)
}

//Multiply by the erf-smoothened coulomb potential in real space, and initialize the index map
inline void multErfCoulomb_thread(int iStart, int iStop, vector3<int> Spad, matrix3<> Rpad,
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
inline void downSample_thread(int iStart, int iStop, const vector3<int>& S, const matrix3<>& GGT,
	const fftw_complex* in, const vector3<int>& Sdense, double* out, double sigma)
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
		out[i] = in[iDense][0] + (4*M_PI) * (Gsq ? (1.-exp(-hlfSigmaSq*Gsq))/Gsq : hlfSigmaSq);
	)
}


CoulombIsolated::CoulombIsolated(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Vc(gInfo)
{
	//Read precomputed kernel from file if supplied
	if(params.filename.length())
	{	FILE* fp = fopen(params.filename.c_str(), "rb");
		if(fp)
		{	matrix3<> R; vector3<int> S; double bw;
			fread(&R, sizeof(matrix3<>), 1, fp);
			fread(&S, sizeof(vector3<int>), 1, fp);
			fread(&bw, sizeof(double), 1, fp);
			if(R != gInfo.R)
				logPrintf("Precomputed coulomb kernel file '%s' has different lattice vectors (recomputing it now)\n", params.filename.c_str());
			else if(!(S == gInfo.S))
				logPrintf("Precomputed coulomb kernel file '%s' has different sample count (recomputing it now)\n", params.filename.c_str());
			else if(bw != params.borderWidth)
				logPrintf("Precomputed coulomb kernel file '%s' has different border width (recomputing it now)\n", params.filename.c_str());
			else if(fread(Vc.data, sizeof(double), gInfo.nG, fp) != unsigned(gInfo.nG))
				logPrintf("Error reading precomputed coulomb kernel from '%s' (computing it now)\n", params.filename.c_str());
			else
			{	logPrintf("Successfully read precomputed coulomb kernel from '%s'\n", params.filename.c_str());
				Vc.set();
				return;
			}
		}
		else logPrintf("Could not open precomputed coulomb kernel file '%s' (computing it now)\n", params.filename.c_str());
	}
	
	//Select gauss-smoothing parameter:
	double maxBorderWidth = 0.5 * ws.inRadius();
	if(params.borderWidth > maxBorderWidth)
		die("Border width %lg bohrs must be less than half the Wigner-Seitz cell in-radius = %lg bohrs.\n",
			params.borderWidth, maxBorderWidth);
	double sigmaBorder = 0.1 * params.borderWidth;
	double sigma = 0.1*ws.inRadius() - sigmaBorder; //so that 10(sigma+sigmaBorder) < inRadius
	logPrintf("Selecting gaussian width %lg bohrs (for border width %lg bohrs).\n", sigmaBorder, params.borderWidth);
	
	//Set up dense integration grids:
	logPrintf("Setting up FFT grids: ");
	double Gnyq = 10./sigmaBorder; //lower bound on Nyquist frequency
	vector3<int> Sdense; //dense fft sample count
	matrix3<> Rpad; vector3<int> Spad; //padded lattice vectors and sample count
	for(int k=0; k<3; k++)
	{	Sdense[k] = std::max(gInfo.S[k], 2*int(ceil(Gnyq * gInfo.R.column(k).length() / (2*M_PI))));
		while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
		//Pad the super cell by border width:
		Spad[k] = Sdense[k] + 2*ceil(Sdense[k]*params.borderWidth/gInfo.R.column(k).length());
		while(!fftSuitable(Spad[k])) Spad[k]+=2;
		Rpad.set_col(k, gInfo.R.column(k) * (Spad[k]*1./Sdense[k]));
	}
	logPrintf("[ %d %d %d ], and [ %d %d %d ] padded.\n",
		Sdense[0], Sdense[1], Sdense[2], Spad[0], Spad[1], Spad[2]);
	int nGdense = Sdense[0] * Sdense[1] * (1+Sdense[2]/2);
	int nGpad = Spad[0] * Spad[1] * (1+Spad[2]/2); //number of symmetry reduced reciprocal lattice vectors
	int nrPad = Spad[0] * Spad[1] * Spad[2]; //number of real space points
	double detRpad = fabs(det(Rpad));
	
	//Construct padded Wigner-Seitz cell, and check border:
	logPrintf("For padded lattice, "); WignerSeitz wsPad(Rpad);
	std::vector<vector3<>> vArr = ws.getVertices();
	matrix3<> invRpad = inv(Rpad);
	for(vector3<> v: vArr)
		if(wsPad.boundaryDistance(wsPad.restrict(invRpad*v)) < 0.9*params.borderWidth)
			die("Padded Wigner-Seitz cell does not fit inside Wigner-Seitz cell of padded lattice.\n"
				"This can happen for reducible lattice vectors; HINT: the reduced lattice vectors,\n"
				"if different from the input lattice vectors, are printed during Symmetry setup.\n");
	
	//Plan Fourier transforms:
	assert(nrPad > 0); //overflow check
	fftw_complex* padArr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nGpad);
	fftw_complex* denseArr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nGdense);
	int* densePadMap = new int[nGpad*2]; //index map from padded to dense grid
	double* padRealArr = (double*)padArr;
	double* denseRealArr = (double*)denseArr;
	if(!padArr || !denseArr || !densePadMap)
		die("Insufficient memory (need %.1fGB). Hint: try increasing border width.\n",
			(nGpad+nGdense)*sizeof(complex)*1e-9 + 2*nGpad*sizeof(int)*1e-9);
	logPrintf("Planning fourier transforms ... "); logFlush();
	fftw_plan_with_nthreads(nProcsAvailable);
	fftw_plan fftPlanC2R = fftw_plan_dft_c2r_3d(Spad[0], Spad[1], Spad[2], padArr, padRealArr, FFTW_ESTIMATE);
	fftw_plan fftPlanR2C = fftw_plan_dft_r2c_3d(Sdense[0], Sdense[1], Sdense[2], denseRealArr, denseArr, FFTW_ESTIMATE);
	logPrintf("Done.\n");
	
	//Initialize smoothed theta function on padded lattice:
	logPrintf("Computing truncation shape function ... "); logFlush();
	std::vector<Simplex<3>> sArr = ws.getSimplices();
	threadLaunch(simplexTheta_thread, nGpad, Spad, (2.0*M_PI)*(~invRpad), detRpad, padArr, &sArr, sigmaBorder);
	fftw_execute(fftPlanC2R);
	logPrintf("Done.\n");
	
	//Multiply by long-ranged erf/r and initialize index map (periodic map via Wigner-Seitz restrictions):
	logPrintf("Applying truncation to Coulomb kernel ... "); logFlush();
	threadLaunch(multErfCoulomb_thread, nrPad, Spad, Rpad, padRealArr, &wsPad, sigma, detRpad/nrPad, densePadMap, Sdense, &ws);
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
	threadLaunch(downSample_thread, gInfo.nG, gInfo.S, gInfo.GGT, denseArr, Sdense, Vc.data, sigma);
	fftw_free(denseArr);
	Vc.set();
	logPrintf("Done.\n");
	
	//Save kernel if requested:
	if(params.filename.length())
	{	logPrintf("Saving isolated coulomb kernel to '%s' ... ", params.filename.c_str()); logFlush();
		FILE* fp = fopen(params.filename.c_str(), "wb");
		if(!fp) die("could not open file for writing.\n");
		fwrite(&gInfo.R, sizeof(matrix3<>), 1, fp);
		fwrite(&gInfo.S, sizeof(vector3<int>), 1, fp);
		fwrite(&params.borderWidth, sizeof(double), 1, fp);
		fwrite(Vc.data, sizeof(double), gInfo.nG, fp);
		fclose(fp);
		logPrintf("Done.\n");
	}
}

DataGptr CoulombIsolated::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombIsolated::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	return EwaldIsolated(gInfo, ws, true, params.borderWidth + params.ionMargin).energyAndGrad(pointCharges);
}


//----------------- class CoulombSpherical ---------------------

CoulombSpherical::CoulombSpherical(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Rc(params.Rc)
{	double RcMax = ws.inRadius();
	if(Rc > RcMax)
		die("Spherical truncation radius %lg exceeds Wigner-Seitz cell in-radius of %lg bohrs.\n", Rc, RcMax);
	if(!Rc) Rc = RcMax;
	logPrintf("Initialized spherical truncation of radius %lg bohrs\n", Rc);
}

DataGptr CoulombSpherical::operator()(DataGptr&& in) const
{	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombSpherical_calc(Rc), in->dataPref(false));
	return in;
}

double CoulombSpherical::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	return EwaldIsolated(gInfo, ws, false, Rc - params.ionMargin).energyAndGrad(pointCharges);
}
