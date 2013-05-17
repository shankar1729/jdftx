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
	diagMatrix width; //exponential widths (this is the 'state' for the minimizer)
	diagMatrix alpha; //atom polarizabilities (inner linear solve at each width set)
	int nAtoms, nSpecies;
	double R0;
	
	SphericalFit(const Everything& e, const ColumnBundle& V, const matrix& X, double alphaTot)
	: e(e), V(V), X(X), alphaTot(alphaTot)
	{
		//Guess initial exponential width based on the vdW R0 parameter:
		logSuspend(); VanDerWaals vdW(e); logResume();
		nAtoms = 0;
		nSpecies = e.iInfo.species.size();
		for(int i=0; i<nSpecies; i++)
		{	const SpeciesInfo& sp = *(e.iInfo.species[i]);
			width.push_back(0.15*vdW.getParams(sp.atomicNumber).R0);
			nAtoms += sp.atpos.size();
		}
		alpha.resize(nSpecies);
		
		//Precompute quantities independent of state:
		matrix VdagKV = V^K(V);
		R0 = 0.5*trace(VdagKV * X * VdagKV * X).real();
	}
	
	void step(const diagMatrix& dir, double alpha)
	{	axpy(alpha, dir, width);
	}
	
	double compute(diagMatrix* grad=0)
	{	//Initialize basis functions for model dielectric matrix at current parameters:
		ColumnBundle U = V.similar(3*nAtoms); //basis functions for the model response
		ColumnBundle dU = U.similar(); //derivative of basis functions w.r.t width parameter
		double dG = 0.02; int nG = int(ceil(e.gInfo.GmaxGrid/dG)) + 5;
		double normFac = sqrt(4*M_PI/3) / e.gInfo.detR;
		std::vector<int> colStart, colStop;
		int colOffset = 0;
		for(int i=0; i<nSpecies; i++)
		{	const SpeciesInfo& sp = *(e.iInfo.species[i]);
			if(width[i] <= 0. || width[i] > 5.) return NAN;
			//Initialize the model radial functions:
			std::vector<double> f_samples(nG), df_samples(nG);
			for(int iG=0; iG<nG; iG++)
			{	double G = iG*dG, wG = width[i]*G;
				f_samples[iG] = normFac * G/pow(1+wG*wG, 3);
				df_samples[iG] = -6.*normFac * G*G*wG/pow(1+wG*wG, 4);
			}
			RadialFunctionG f, df;
			f.init(1, f_samples, dG);
			df.init(1, df_samples, dG);
			//Set the basis functions:
			const Basis& basis = *(U.basis);
			int l = 1;
			int atomStride = basis.nbasis;
			colStart.push_back(colOffset);
			for(int m=-l; m<=l; m++)
			{	callPref(Vnl)(basis.nbasis, atomStride, sp.atpos.size(), l, m, vector3<>(), basis.iGarrPref, e.gInfo.G, sp.atposPref, f, U.dataPref()+colOffset*basis.nbasis, false, vector3<complex*>());
				callPref(Vnl)(basis.nbasis, atomStride, sp.atpos.size(), l, m, vector3<>(), basis.iGarrPref, e.gInfo.G, sp.atposPref, df, dU.dataPref()+colOffset*basis.nbasis, false, vector3<complex*>());
				colOffset += sp.atpos.size();
			}
			colStop.push_back(colOffset);
			f.free();
			df.free();
		}
		
		matrix VdagKU = V^K(U), UdagKV = dagger(VdagKU), UdagKU = U^K(U);
		diagMatrix diag_UdagKV_X_VdagKU = diag(UdagKV * X * VdagKU);
		
		//Solve for the atom polarizabilities at current width: (constrained linear solve)
		matrix A(nSpecies,nSpecies), b(nSpecies,1), N(nSpecies,1);
		for(int i=0; i<nSpecies; i++)
		{	b.data()[b.index(i,0)] = trace(diag_UdagKV_X_VdagKU(colStart[i],colStop[i]));
			N.data()[N.index(i,0)] = (1./3) * (colStop[i] - colStart[i]);
			for(int j=0; j<nSpecies; j++)
			{	matrix UdagKUsub = UdagKU(colStart[i],colStop[i], colStart[j],colStop[j]);
				A.data()[A.index(i,j)] = (1./3) * trace(UdagKUsub * dagger(UdagKUsub)).real();
			}
		}
		matrix invA = inv(A);
		double lambda = ((trace(dagger(N)*invA*b) + alphaTot) / trace(dagger(N)*invA*N)).real(); //Lagrange multiplier for net polarizability constraint
		matrix alphaMat = invA * (lambda*N - b);
		//Convert above solution to diagonal matrix in the model basis (and store to this->alpha):
		diagMatrix Y(3*nAtoms);
		for(int i=0; i<nSpecies; i++)
		{	alpha[i] = alphaMat.data()[alphaMat.index(i,0)].real();
			std::fill(Y.data()+colStart[i], Y.data()+colStop[i], (1./3)*alpha[i]);
		}
		
		//Compute residual and gradient:
		double R = R0 + 0.5*trace(UdagKU * Y * UdagKU * Y).real() + trace(diag_UdagKV_X_VdagKU * Y);
		if(grad)
		{	diagMatrix R_w = diag(K(dU) ^ (V * (X * VdagKU * Y) + U * (Y * UdagKU * Y)));
			grad->resize(nSpecies);
			for(int i=0; i<nSpecies; i++)
				grad->at(i) = 2. * trace(R_w(colStart[i],colStop[i]));
			constrain(*grad);
		}
		return R;
	}

	static inline void K_thread(int bStart, int bStop, const Everything* e, vector3<> dk, const ColumnBundle* rho, ColumnBundle* Krho)
	{	for(int b=bStart; b<bStop; b++)
			Krho->setColumn(b, e->gInfo.detR * (*(e->coulomb))(rho->getColumn(b), dk, 0.));
	}
	ColumnBundle K(ColumnBundle Y)
	{	ColumnBundle KY = Y.similar();
		suspendOperatorThreading();
		threadLaunch(K_thread, Y.nCols(), &e, vector3<>(), &Y, &KY);
		resumeOperatorThreading();
		return KY;
	}
	
	void print()
	{	logPrintf("\tSpherical polarizability parameters:\n");
		for(unsigned i=0; i<e.iInfo.species.size(); i++)
			logPrintf("\t   Site '%s'   alpha: %lf   width: %lf\n",
				e.iInfo.species[i]->name.c_str(), alpha[i], width[i]);
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
	mp.nDim = e.iInfo.species.size();
	mp.linePrefix = "\tChiSphericalFit: ";
	mp.energyLabel = "residual";
	mp.nIterations = 100;
	mp.energyDiffThreshold = 1e-8;
	SphericalFit fit(e, V, Xext, alpha);
	fit.minimize(mp);
	fit.print();
	
	finalizeSystem();
	return 0;
}
