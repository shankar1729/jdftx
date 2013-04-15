/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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
#include <electronic/Polarizability.h>
#include <electronic/VanDerWaals.h>
#include <electronic/operators.h>
#include <electronic/SpeciesInfo_internal.h>
#include <commands/parser.h>
#include <core/DataMultiplet.h>
#include <core/Util.cpp>


struct SphericalFit : public Minimizable<diagMatrix>
{
	const Everything& e;
	const ColumnBundle& V;
	const matrix& X;
	double alphaTot;
	diagMatrix state;
	int nAtoms;
	double R0;
	
	SphericalFit(const Everything& e, const ColumnBundle& V, const matrix& X, double alphaTot)
	: e(e), V(V), X(X), alphaTot(alphaTot)
	{
		//Guess initial parameters based on the vdW C6 / R0 parameters:
		logSuspend(); VanDerWaals vdW(e); logResume();
		double C6tot = 0.;
		nAtoms = 0;
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
		{	const SpeciesInfo& sp = *(e.iInfo.species[i]);
			VanDerWaals::AtomParams params = vdW.getParams(sp.atomicNumber);
			state.push_back(0.15*params.R0); //this is the typical density decay length for atoms
			state.push_back(params.C6); //guess polarizabilities in proportion to C6
			C6tot += params.C6 * sp.atpos.size();
			nAtoms += sp.atpos.size();
		}
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
			state[2*i+1] *= (alphaTot / C6tot);
		
		//Precompute quantities independent of state:
		matrix VdagOV = V^O(V);
		R0 = 0.5*trace(VdagOV * X * VdagOV * X).real();
	}
	
	void step(const diagMatrix& dir, double alpha)
	{	axpy(alpha, dir, state);
	}
	
	double compute(diagMatrix* grad=0)
	{	//Initialize model dielectric matrix at current parameters:
		diagMatrix Y(3*nAtoms); //Diagonal entries for the model response
		ColumnBundle U = V.similar(3*nAtoms); //basis functions for the model response
		ColumnBundle dU = U.similar(); //derivative of basis functions w.r.t width parameter
		double dG = 0.02; int nG = int(ceil(e.iInfo.GmaxLoc/dG)) + 5;
		double normFac = sqrt(4*M_PI/3) / e.gInfo.detR;
		int colOffset = 0;
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
		{	const SpeciesInfo& sp = *(e.iInfo.species[i]);
			const double& width = state[2*i];
			const double& alpha = state[2*i+1];
			if(width <= 0. || width > 5.) return NAN;
			//Set the magnitudes:
			std::fill_n(Y.begin()+colOffset, 3*sp.atpos.size(), (1./3)*alpha);
			//Initialize the model radial functions:
			std::vector<double> f_samples(nG), df_samples(nG);
			for(int i=0; i<nG; i++)
			{	double G = i*dG, wG = width*G;
				f_samples[i] = normFac * G/pow(1+wG*wG, 3);
				df_samples[i] = -6.*normFac * G*G*wG/pow(1+wG*wG, 4);
			}
			RadialFunctionG f, df;
			f.init(1, f_samples, dG);
			df.init(1, df_samples, dG);
			//Set the basis functions:
			const Basis& basis = *(U.basis);
			int l = 1;
			int atomStride = basis.nbasis;
			for(int m=-l; m<=l; m++)
			{	callPref(Vnl)(basis.nbasis, atomStride, sp.atpos.size(), l, m, vector3<>(), basis.iGarrPref, e.gInfo.G, sp.atposPref, f, U.dataPref()+colOffset*basis.nbasis, false, vector3<complex*>());
				callPref(Vnl)(basis.nbasis, atomStride, sp.atpos.size(), l, m, vector3<>(), basis.iGarrPref, e.gInfo.G, sp.atposPref, df, dU.dataPref()+colOffset*basis.nbasis, false, vector3<complex*>());
				colOffset += sp.atpos.size();
			}
			f.free();
		}
		
		//Compute residual and gradient:
		matrix VdagOU = V^O(U), UdagOV = dagger(VdagOU), UdagOU = U^O(U);
		double R = R0 + 0.5*trace(UdagOU * Y * UdagOU * Y).real() + trace(UdagOV * X * VdagOU * Y).real();
		if(grad)
		{	diagMatrix R_Y = diag(UdagOV * X * VdagOU + UdagOU * Y * UdagOU);
			diagMatrix R_w = diag(O(dU) ^ (V * (X * VdagOU * Y) + U * (Y * UdagOU * Y)));
			int colStart = 0;
			grad->resize(state.size());
			for(unsigned i=0; i<e.iInfo.species.size(); i++)
			{	int colStop = colStart + e.iInfo.species[i]->atpos.size()*3;
				grad->at(2*i) = 2. * trace(R_w(colStart,colStop));
				grad->at(2*i+1) = (1./3) * trace(R_Y(colStart,colStop));
				colStart = colStop;
			}
			constrain(*grad);
		}
		return R;
	}

	//project out changes to total alpha:
	void constrain(diagMatrix& dir)
	{	double dalphaMean = 0.;
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
			dalphaMean += e.iInfo.species[i]->atpos.size() * dir[2*i+1];
		dalphaMean /= nAtoms;
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
			dir[2*i+1] -= dalphaMean;
	}

	void print()
	{	logPrintf("\tSpherical polarizability parameters:\n");
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
		{	const SpeciesInfo& sp = *(e.iInfo.species[i]);
			const double& width = state[2*i];
			const double& alpha = state[2*i+1];
			logPrintf("\t   Site '%s'   alpha: %lf   width: %lf\n", sp.name.c_str(), alpha, width);
		}
	}
};


inline void rArrSet(size_t iStart, size_t iStop, vector3<int> S, matrix3<> R, vector3<double*> rArr)
{	matrix3<> h; for(int k=0; k<3; k++) h.set_col(k, (1./S[k])*R.column(k)); //mesh vectors
	THREAD_fullGspaceLoop
	(	storeVector(h * iG, rArr, i);
	)
}

int main(int argc, char** argv)
{	if(argc != 2) { logPrintf("Usage: SphericalChi <jdftx-input-file>\n"); return 1; }
	initSystem(argc, argv);
	
	//Read the input file:
	Everything e;
	parse(argv[1], e);
	//Skip reading the wavefunctions and initialzing empty states to save time:
	assert(e.dump.polarizability);
	Polarizability& pol = *e.dump.polarizability;
	int nV = e.eInfo.nElectrons/2;
	int nC = e.eInfo.nBands - nV;
	if(pol.Ecut<=0.) pol.Ecut = 4.*e.cntrl.Ecut;
	if(pol.nEigs<=0) pol.nEigs = nV * nC;
	e.eInfo.nBands = nV;
	e.eVars.wfnsFilename.clear();
	//Perform JDFTx initialization:
	e.setup();
	
	logPrintf("----- Processing polarizability eigenfunctions -----\n");
	
	//Setup the basis and read in polarizability:
	logPrintf("Setting up reduced basis at Ecut=%lg: ", pol.Ecut);
	Basis basis; basis.setup(e.gInfo, e.iInfo, pol.Ecut, pol.dk);
	logPrintf("Reading polarizability eigensystem\n");
	ColumnBundle V(pol.nEigs, basis.nbasis, &basis);
	matrix Xext(pol.nEigs, pol.nEigs);
	((ManagedMemory&)V).read(e.dump.getFilename("pol_basis").c_str());
	Xext.read(e.dump.getFilename("pol_Xext").c_str());
	
	logPrintf("Computing eigenfunction dipole moments\n");
	matrix P;
	{	//Create array of r in real space:
		DataRptrVec rArr(e.gInfo);
		threadLaunch(rArrSet, e.gInfo.nr, e.gInfo.S, e.gInfo.R, rArr.data());
		
		//Convert to a columnbundle projector:
		ColumnBundle OJr(3, basis.nbasis, &basis);
		for(int k=0; k<3; k++) OJr.setColumn(k, O(J(Complex(rArr[k]))));
	
		P = OJr ^ V;
		string fname = e.dump.getFilename("pol_P");
		logPrintf("\tDumping %s ... ", fname.c_str()); logFlush();
		P.write(fname.c_str());
		logPrintf("Done!\n"); logFlush();
	}
	double alpha = (-1./3) * trace(P * Xext * dagger(P)).real(); //isotropic polarizability
	logPrintf("\tIsotropic polarizability: %lg bohr^3\n", alpha);
	
	
	logPrintf("Fitting polarizability to spherical atom-centered model\n");
	MinimizeParams mp;
	mp.nDim = 2. * e.iInfo.species.size();
	mp.linePrefix = "\tChiSphericalFit: ";
	mp.energyLabel = "Residual";
	mp.nIterations = 100;
	mp.knormThreshold = 1e-6;
	SphericalFit fit(e, V, Xext, alpha);
	for(int i=0; i<10; i++) fit.minimize(mp);
	fit.print();
	
	finalizeSystem();
	return 0;
}
