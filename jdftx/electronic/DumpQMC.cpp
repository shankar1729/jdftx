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

#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/ElecMinimizer.h>
#include <core/Blip.h>
#include <core/matrix.h>
#include <core/Thread.h>
#include <core/Operators.h>
#include <core/ScalarFieldIO.h>
#include <config.h>
#include <map>

int nAtomsTot(const IonInfo& iInfo)
{	unsigned res=0;
	for(auto sp: iInfo.species) res += sp->atpos.size();
	return res;
}

struct ImagPartMinimizer: public Minimizable<ElecGradient>  //Uses only the Haux entries of ElecGradient
{	const Everything& e;
	std::vector<matrix> mask; //matrices of ones and zeroes indicating which columns are allowed to mix
	std::vector<matrix> U, Haux; //unitary rotations and the hermitian matrices that generate them
	int nDim; //dimension of minimize (no minimize needed if zero)
	
	ImagPartMinimizer(const Everything& e) : e(e), mask(e.eInfo.nStates), U(e.eInfo.nStates), Haux(e.eInfo.nStates)
	{	//Initialize mask:
		nDim = 0;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	Haux[q] = zeroes(e.eInfo.nBands, e.eInfo.nBands);
			mask[q] = zeroes(e.eInfo.nBands, e.eInfo.nBands);
			const diagMatrix& eigs = e.eVars.Hsub_eigs[q];
			complex* maskData = mask[q].data();
			for(int b1=0; b1<e.eInfo.nBands; b1++)
				for(int b2=0; b2<e.eInfo.nBands; b2++)
					if(fabs(eigs[b1]-eigs[b2]) < 1e-4)
					{	maskData[mask[q].index(b1,b2)] = 1.;
						if(b1<b2) nDim++;
					}
		}
		mpiWorld->allReduce(nDim, MPIUtil::ReduceSum);
	}
	
	void step(const ElecGradient& dir, double alpha)
	{	assert(dir.Haux.size() == Haux.size());
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			if(dir.Haux[q]) axpy(alpha, dir.Haux[q], Haux[q]);
	}
	
	double compute(ElecGradient* grad, ElecGradient* Kgrad)
	{	if(grad) grad->init(e);
		double imagErr = 0.;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	matrix Bevecs; diagMatrix Beigs;
			U[q] = cis(Haux[q], &Bevecs, &Beigs);
			ColumnBundle Crot = e.eVars.C[q] * U[q], imagErr_Crot;
			if(grad) imagErr_Crot = Crot.similar();
			for(int b=0; b<e.eInfo.nBands; b++)
			{	ScalarField imagIpsi = Imag(I(Crot.getColumn(b,0)));
				imagErr += e.gInfo.dV * dot(imagIpsi,imagIpsi);
				if(grad)
					imagErr_Crot.setColumn(b,0, Idag(Complex(0*imagIpsi, 2*e.gInfo.dV*imagIpsi)));
			}
			if(grad)
				grad->Haux[q] = dagger_symmetrize(cis_grad(U[q] * (imagErr_Crot ^ Crot) * dagger(U[q]), Bevecs, Beigs));
		}
		mpiWorld->allReduce(imagErr, MPIUtil::ReduceSum);
		if(grad)
		{	constrain(*grad);
			if(Kgrad) *Kgrad = *grad;
		}
		return imagErr;
	}
	
	void constrain(ElecGradient& grad)
	{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			callPref(eblas_zmul)(mask[q].nData(), mask[q].dataPref(),1, grad.Haux[q].dataPref(),1); //apply mask
	}
	
	double sync(double x) const { mpiWorld->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
};

void Dump::dumpQMC()
{
	const IonInfo &iInfo = e->iInfo;
	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	const Energies &ener = e->ener;
	const Basis& basis = e->basis[0];
	const GridInfo& gInfo = e->gInfo;

	if(e->eInfo.isNoncollinear())
		die("QMC output is not supported for noncollinear spin modes.\n");
	
	BlipConverter blipConvert(gInfo.S);
	string fname; ofstream ofs;
	
	ScalarFieldTilde nTilde = J(eVars.get_nTot());
	
	//Total effective potential on electrons due to fluid:
	ScalarField Vdiel = eVars.fluidParams.fluidType==FluidNone ? 0 : I(eVars.d_fluid + eVars.V_cavity);
	nullToZero(Vdiel, gInfo);
	
	//QMC energy corrections (all information needed for different self-consistency schemes):
	const double& A_diel = e->ener.E["A_diel"];
	double Atilde_diel = A_diel - dot(nTilde,O(J(Vdiel))); //Legendre transform to get final correction
	double Anuc_diel = eVars.d_fluid ? dot(eVars.d_fluid, O(iInfo.rhoIon)) : 0; //interaction with nuclear charge
	logPrintf("QMC energy correction:\n");
	logPrintf("\tA_diel      = %25.16lf\n", A_diel);
	logPrintf("\tAtilde_diel = %25.16lf\n", Atilde_diel);
	logPrintf("\tAnuc_diel   = %25.16lf\n", Anuc_diel);

	//-------------------------------------------------------------------------------------------
	//Output potential directly in BLIP-function basis (cubic B-splines)
	ScalarField VdielBlip(blipConvert(Vdiel));

	#define StartDump(varName) \
		fname = getFilename(varName); \
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush(); \
		if(!mpiWorld->isHead()) fname = "/dev/null";
	StartDump("expot.data")
	ofs.open(fname);
	ofs.precision(12);
	ofs.setf(std::ios::scientific);
	ofs <<
		"START HEADER\n"
		" CASINO Blip external potential exported by " PACKAGE_NAME " " VERSION_MAJOR_MINOR_PATCH "\n"
		"END HEADER\n"
		"\n"
		"START VERSION\n"
		" 2\n"
		"END VERSION\n"
		"\n"
		"START EXPOT\n"
		" Title\n"
		"  Potential due to fluid in JDFT\n"
		" Number of sets\n"
		"  1\n"
		" START SET 1\n"
		"  Particle types affected\n"
		"   All\n"
		"  Periodicity\n"
		"   PERIODIC\n"
		"  Type of representation\n"
		"   BLIP_PERIODIC\n"
		"  Direction\n"
		"   X\n"
		"  Number of such potentials to add\n"
		"   1\n"
		"  Origin for potential\n"
		"   0.0 0.0 0.0\n"
		"  Potential\n"
		"   Blip grid size\n"
		"    " << gInfo.S[0] << " " << gInfo.S[1] << " " << gInfo.S[2] << "\n"
		"   Blip coefficients\n";
	double* VdielBlipData = VdielBlip->data();
	for(int i=0; i<gInfo.nr; i++) ofs << VdielBlipData[i] << "\n";
	ofs <<
		" END SET 1\n"
		"END EXPOT\n";
	ofs.close();
	logPrintf("done.\n"); logFlush();

	//Output Vexternal in BLIPs (if present)
	for(unsigned s=0; s<eVars.Vexternal.size(); s++)
	{	string varName = "VexternalBlip";
		if(eVars.n.size()==2)
			varName += s==0 ? "Up" : "Dn";
		fname = getFilename(varName);
		logPrintf("Dumping '%s'...", fname.c_str()); logFlush();
		if(mpiWorld->isHead()) saveRawBinary(blipConvert(eVars.Vexternal[s]), fname.c_str());
		logPrintf("done.\n"); logFlush();
	}
	
	//-------------------------------------------------------------------------------------------
	//Output wavefunctions directly in BLIP-function basis (cubic B-splines)
	StartDump("bwfn.data")
	logPrintf("\n");
	ofs.open(fname);
	ofs.precision(12);
	ofs.setf(std::ios::scientific);
	ofs <<
		"CASINO Blip orbitals exported by  " PACKAGE_NAME "\n"
		"\n"
		"BASIC INFO\n"
		"----------\n"
		" Generated by:\n"
		"  " PACKAGE_NAME " " VERSION_MAJOR_MINOR_PATCH "\n"
		" Method:\n"
		"  DFT\n"
		" DFT Functional:\n"
		"  " << e->exCorr.getName() << "\n"
		" Pseudopotential\n"
		"  Norm-conserving Kleinman-Bylander\n"
		" Plane wave cutoff (au)\n"
		"  " << e->cntrl.Ecut << "\n"
		" Spin polarized:\n"
		"  " << (eInfo.spinType==SpinNone ? 'F' : 'T') << "\n"
		" Total energy (au per primitive cell)\n"
		"  " << ener.E["Etot"] << "\n"
		" Kinetic energy (au per primitive cell)\n"
		"  " << ener.E["KE"] << "\n"
		" Local potential energy (au per primitive cell)\n"
		"  " << ener.E["Eloc"] + A_diel - Atilde_diel << "\n"
		" Non local potential energy (au per primitive cell)\n"
		"  " << ener.E["Enl"] + ener.E["Exc"] << "\n"
		" Electron electron energy (au per primitive cell)\n"
		"  " << ener.E["EH"] << "\n"
		" Ion-ion energy (au per primitive cell)\n"
		"  " << ener.E["Eewald"] + Atilde_diel << "\n"
		" Number of electrons per primitive cell\n"
		"  " << int(eInfo.nElectrons+0.5) << "\n"
		"\n"
		"GEOMETRY\n"
		"--------\n"
		" Number of atoms per primitive cell\n"
		"  " << nAtomsTot(e->iInfo) << "\n"
		" Atomic number and position of the atoms(au)\n";
	for(auto sp: iInfo.species)
		for(unsigned a=0; a<sp->atpos.size(); a++)
		{
			assert(sp->atomicNumber);
			vector3<> coord(gInfo.R * sp->atpos[a]);
			ofs << "  " << sp->atomicNumber << " "
				<< coord[0] << " " << coord[1] << " " << coord[2]
				<< "\n";
		}
	ofs <<
		" Primitive lattice vectors (au)\n"
		"  " << gInfo.R(0,0) << " " << gInfo.R(1,0) << " " << gInfo.R(2,0) << "\n"
		"  " << gInfo.R(0,1) << " " << gInfo.R(1,1) << " " << gInfo.R(2,1) << "\n"
		"  " << gInfo.R(0,2) << " " << gInfo.R(1,2) << " " << gInfo.R(2,2) << "\n"
		"\n"
		"G VECTORS\n"
		"---------\n"
		" Number of G-vectors\n"
		"  " << basis.nbasis << "\n"
		" Gx Gy Gz (au)\n";
	vector3<> Gvec;
	ofs << "  " << Gvec[0] << " " << Gvec[1] << " " << Gvec[2] << "\n"; //G=0
	for(const vector3<int>& iG: basis.iGarr) //G!=0:
		if(iG.length_squared())
		{	vector3<> Gvec = iG * gInfo.G;		
			ofs << "  " << Gvec[0] << " " << Gvec[1] << " " << Gvec[2] << "\n";
		}
	ofs <<
		" Blip Grid\n"
		"  " << gInfo.S[0] << " " << gInfo.S[1] << " " << gInfo.S[2] << "\n"
		"\n"
		"WAVE FUNCTION\n"
		"-------------\n";
	int nSpins = eInfo.nSpins();
	int nkPoints = eInfo.nStates / nSpins;
	ofs <<
		" Number of k-points\n"
		"  " << nkPoints << "\n";

	//Handle degeneracy in Gamma-point only calculations (necessary for removePhase to work correctly)
	std::vector<matrix> Udeg; //rotations to make wavefunctions real
	if(nkPoints==1)
	{	ImagPartMinimizer imin(*e);
		if(imin.nDim)
		{	logPrintf("\tFinding linear combinations of orbitals (degrees of freedom: %d) that are real:\n", imin.nDim);
			MinimizeParams mp;
			mp.nDim = imin.nDim;
			mp.energyLabel = "ImagNorm";
			mp.linePrefix = "\t\tImagMinimize: ";
			mp.energyDiffThreshold = 1e-12;
			mp.dirUpdateScheme = MinimizeParams::LBFGS;
			mp.nIterations = 200;
			mp.fpLog = globalLog;
			imin.minimize(mp);
			Udeg = imin.U;
		}
	}
	
	for(int ik=0; ik<nkPoints; ik++)
	{
		vector3<> k(eInfo.qnums[ik].k * gInfo.G); //cartesian k-vector
		ofs <<
			" k-point # ; # of bands (up spin/down spin) ; k-point coords (au)\n"
			"  " << ik+1 << " " << eInfo.nBands << " " << (nSpins==2 ? eInfo.nBands : 0)
			<< " " << k[0] << " " << k[1] << " " << k[2] << "\n";
		for(int s=0; s<nSpins; s++)
		{
			int spinIndex = 1-2*s;
			int q = ik + nkPoints*s; //net quantum number
			
			//Get relevant wavefunctions and eigenvalues (from another process if necessary)
			ColumnBundle CqTemp; diagMatrix Hsub_eigsqTemp;
			const ColumnBundle* Cq=0; const diagMatrix* Hsub_eigsq=0;
			if(mpiWorld->isHead())
			{	if(eInfo.isMine(q))
				{	Cq = &eVars.C[q];
					Hsub_eigsq = &eVars.Hsub_eigs[q];
				}
				else
				{	Cq = &CqTemp;
					CqTemp.init(eInfo.nBands, e->basis[q].nbasis, &e->basis[q], &eInfo.qnums[q]);
					mpiWorld->recvData(CqTemp, eInfo.whose(q), q);
					Hsub_eigsq = &Hsub_eigsqTemp;
					Hsub_eigsqTemp.resize(eInfo.nBands);
					mpiWorld->recvData(Hsub_eigsqTemp, eInfo.whose(q), q);
					if(Udeg.size())
					{	Udeg[q].init(eInfo.nBands, eInfo.nBands);
						mpiWorld->recvData(Udeg[q], eInfo.whose(q), q);
					}
				}
			}
			else
			{	if(eInfo.isMine(q))
				{	mpiWorld->sendData(eVars.C[q], 0, q);
					mpiWorld->sendData(eVars.Hsub_eigs[q], 0, q);
					if(Udeg.size()) mpiWorld->sendData(Udeg[q], 0, q);
				}
				continue; //only head performs computation below
			}
			
			//Apply degeneracy rotations (if any):
			if(Udeg.size())
			{	CqTemp = (*Cq) * Udeg[q];
				Cq = &CqTemp;
			}
			
			//Loop over bands
			for (int b=0; b < eInfo.nBands; b++)
			{
				logPrintf("\tProcessing state %3d band %3d:\n", q, b); logFlush();
				int bandIndex = b+1 + s*eInfo.nBands;
				ofs <<
					" Band, spin, eigenvalue (au), localized\n"
					"  " << bandIndex << " " << spinIndex << " " << Hsub_eigsq->at(b) << " F\n"
					" " << (nkPoints==1 ? "Real" : "Complex") << " blip coefficients for extended orbitals\n";
				//Get orbital in real space
				complexScalarField phi = I(Cq->getColumn(b,0));
				//Compute kinetic and potential energy (in Vdiel) of original PW orbitals:
				double Tpw = -0.5*dot(phi, Jdag(L(J(phi)))).real();
				double Vpw = gInfo.dV*dot(phi, Vdiel*phi).real();
				//Adjust phase, convert to blip and output in appropriate format:
				if(nkPoints==1)
				{	//Convert orbital to real:
					double phaseMean, phaseSigma, imagErrorRMS;
					removePhase(gInfo.nr, phi->data(), phaseMean, phaseSigma, imagErrorRMS);
					logPrintf("\t\tPhase = %lf +/- %lf\n", phaseMean, phaseSigma); logFlush();
					logPrintf("\t\tImagErrorRMS = %le\n", imagErrorRMS); logFlush();
					//Convert PW to blip and output:
					phi = blipConvert(phi);
					complex* phiData = phi->data();
					for(int i=0; i<gInfo.nr; i++)
						ofs << "  " << phiData[i].real() << "\n";
				}
				else
				{	//Convert PW to blip:
					phi = blipConvert(phi); complex* phiData = phi->data();
					//Output:
					for(int i=0; i<gInfo.nr; i++)
						ofs << "  (" << phiData[i].real() << "," << phiData[i].imag() << ")\n";
				}
				//Compare PW and Blip kinetic and potential energies
				double tMax; int i0max, i1max, i2max;
				double Tblip = ::Tblip(phi, &tMax, &i0max, &i1max, &i2max);
				logPrintf("\t\tKinetic Energy    (PW)    = %.12le\n", Tpw);
				logPrintf("\t\tKinetic Energy    (Blip)  = %.12le (Ratio = %.6lf)\n", Tblip, Tblip/Tpw);
				logPrintf("\t\tMax local KE      (Blip)  = %.12le at cell#(%d,%d,%d)\n", tMax, i0max, i1max, i2max);
				logFlush();
				if(A_diel)
				{	//Vdiel potential contribution only when Adiel is non-zero
					double Vblip = ::Vblip(phi, VdielBlip);
					logPrintf("\t\tInt Vdiel.|phi^2| (PW)    = %.12le\n", Vpw);
					logPrintf("\t\tInt Vdiel.|phi^2| (Blip)  = %.12le (Ratio = %.6lf)\n", Vblip, Vblip/Vpw);
					logFlush();
				}
			}
		}
	}
	logPrintf("\tDone.\n"); logFlush();
	ofs.close();
}
