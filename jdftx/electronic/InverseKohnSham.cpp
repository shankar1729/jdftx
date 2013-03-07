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

#include <electronic/InverseKohnSham.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <electronic/operators.h>
#include <electronic/matrix.h>

off_t fileSize(const char *filename); //defined in Util.cpp

//Compute a guess for inv(chi) using wavefunctions, fillings and eigenvalues of a similar system
class InvertChi : public LinearSolvable<DataRptrCollection>
{
public:
	InvertChi(const Everything& e);
	DataRptrCollection hessian(const DataRptrCollection&);
	DataRptrCollection precondition(const DataRptrCollection& v);
private:
	const Everything& e;
	std::vector<ColumnBundle> C;
	std::vector<diagMatrix> F;
	std::vector<diagMatrix> eigs;
};


InverseKohnSham::InverseKohnSham(Everything& e) : e(e), n(e.eVars.n.size()), gaussCutoff(e.gInfo), invertChi(0)
{	//Setup minimization parameters:
	e.inverseKSminParams.nDim = e.gInfo.nr * n.size();
	e.inverseKSminParams.fpLog = globalLog;
	e.inverseKSminParams.linePrefix = "InverseKohnShamMinimize: ";
	e.inverseKSminParams.energyLabel = "minusWs";
	
	if(!e.cntrl.invertKS_nonlocal) //Eliminate non-local pseudopotentials:
		e.iInfo.species.clear(); //ions have served their purpose (computing nElectrons, symmetries and Vlocps)
	
	//Prepare initial state (if invalid):
	if(!e.eVars.Vexternal.size())
		nullToZero(e.eVars.Vexternal, e.gInfo, n.size());
	
	initGaussianKernel(gaussCutoff, e.cntrl.invertKS_sigma); //cutoff kernel
	
	//Initialize inv(chi) preconditioner if specified:
	if(e.cntrl.invertKS_chiGuessFilename.length())
		invertChi = std::make_shared<InvertChi>(e);
}

void InverseKohnSham::step(const DataRptrCollection& dir, double alpha)
{	axpy(alpha, dir, e.eVars.Vexternal);
	e.eVars.EdensityAndVscloc(e.ener); //update Vscloc for subsequent band structure minimize
}

double InverseKohnSham::compute(DataRptrCollection* grad)
{	logPrintf("\n----------- Band structure minimization -------------\n"); logFlush();
	elecMinimize(e);
	//---Compute inverse Kohn sham functional, Ws [J Chem Phys 118, 2498 (2003)]
	//Compute the density from the band-structure solve:
	for(unsigned s=0; s<n.size(); s++) initZero(n[s], e.gInfo); //Initialize to zero
	for(int q=0; q<e.eInfo.nStates; q++)
		n[e.eInfo.qnums[q].index()] += e.eInfo.qnums[q].weight * diagouterI(e.eVars.F[q], e.eVars.C[q]);
	for(unsigned s=0; s<n.size(); s++) e.symm.symmetrize(n[s]);
	//Compute the energy of the eigenfunctions in the net local potential:
	double minusWs = dot(e.eVars.Vscloc, e.eVars.n);
	for(int q=0; q<e.eInfo.nStates; q++)
		minusWs -= e.eInfo.qnums[q].weight * trace(e.eVars.F[q] * e.eVars.Hsub_eigs[q]);
	if(grad)
	{	grad->resize(n.size());
		for(unsigned s=0; s<n.size(); s++)
		{	(*grad)[s] = e.gInfo.dV * (e.eVars.n[s] - n[s]);
			e.symm.symmetrize((*grad)[s]);
		}
	}
	return minusWs;
}


DataRptrCollection InverseKohnSham::precondition(const DataRptrCollection& grad)
{	DataRptrCollection Kgrad = clone(grad);
	//Right cutoff:
	constrain(Kgrad);
	for(unsigned s=0; s<n.size(); s++) Kgrad[s] = I(gaussCutoff*J(Kgrad[s]));
	
	if(invertChi) //precondition using chi of some (hopefully similar) electronic system
	{	MinimizeParams invertChiMinParams = e.inverseKSminParams;
		invertChiMinParams.linePrefix = "InverseChiLinearCG: ";
		invertChiMinParams.knormThreshold = 0;
		invertChiMinParams.nIterations = 20;
		nullToZero(invertChi->state, e.gInfo, grad.size());
		invertChi->solve(Kgrad, invertChiMinParams);
		Kgrad = invertChi->state;
		invertChi->state.clear();
	}
	//Left cutoff:
	for(unsigned s=0; s<n.size(); s++) Kgrad[s] = I(gaussCutoff*J(Kgrad[s]));
	constrain(Kgrad);
	return Kgrad;
}

void InverseKohnSham::constrain(DataRptrCollection& dir)
{	for(unsigned s=0; s<dir.size(); s++)
	{	e.symm.symmetrize(dir[s]); //symmetrize
		dir[s] += (-1./e.gInfo.nr) * sum(dir[s]); //project out G=0
	}
}

bool InverseKohnSham::report(int iter)
{	e.dump(DumpFreq_Ionic, iter); //Each step is an electronic solve, so this roughly equivalent to an ionic step in cost
	return false;
}

//----------------- Implementation of class InvertChi ------------------

InvertChi::InvertChi(const Everything& e) : e(e)
{
	string fname;
	#define setFname(varName) \
		fname = e.cntrl.invertKS_chiGuessFilename; \
		fname.replace(fname.find("$VAR"),4, #varName);
	
	//Read eigenvalues:
	setFname(eigenvals);
	FILE* fp = fopen(fname.c_str(), "r");
	if(!fp) die("Could not open %s for reading chiGuess eigenvalues.\n", fname.c_str());
	std::vector<double> eigArr;
	while(!feof(fp))
	{	double eig;
		if(fscanf(fp, "%lf", &eig)==1)
			eigArr.push_back(eig);
	}
	fclose(fp);
	if(eigArr.size() % e.eInfo.nStates != 0)
		die("Number of eigenvalues in %s = %lu is not a multiple of nStates = %d.\n",
			fname.c_str(), eigArr.size(), e.eInfo.nStates);
	int nBandsChi = eigArr.size() / e.eInfo.nStates;
	eigs.resize(e.eInfo.nStates);
	for(int q=0; q<e.eInfo.nStates; q++)
		eigs[q].assign(eigArr.begin()+q*nBandsChi, eigArr.begin()+(q+1)*nBandsChi);
	
	//Read / compute fillings:
	setFname(fillings)
	fp = fopen(fname.c_str(), "r");
	F.resize(e.eInfo.nStates);
	if(fp) //read from file
	{	for(int q=0; q<e.eInfo.nStates; q++)
		{	F[q].resize(nBandsChi);
			F[q].scan(fp);
			F[q] *= (0.5*e.eVars.n.size()); //Internal fillings 0 to 1
		}
		fclose(fp);
	}
	else //use eVars.F
	{	for(int q=0; q<e.eInfo.nStates; q++)
		{	F[q] = e.eVars.F[q];
			F[q].resize(nBandsChi, 0.); //pad with zeros (or drop extra entries)
		}
	}
	
	//Read wavefunctions:
	setFname(wfns)
	init(C, e.eInfo.nStates, nBandsChi, &e.basis[0], &e.eInfo.qnums[0]);
	size_t expectedLen = 0;
	for(const ColumnBundle& psi: C)
		expectedLen += psi.nData() * sizeof(complex);
	size_t fileLen = fileSize(fname.c_str());
	if(fileLen<0)
		die("Error accessing chiGuess wavefunctions %s.\n", fname.c_str());
	if(fileLen!=expectedLen)
		die("Length of chiGuess wavefunction file %s is %lu (expected %lu bytes).\n",
			fname.c_str(), fileLen, expectedLen);
	fp = fopen(fname.c_str(), "rb");
	for(int q=0; q<e.eInfo.nStates; q++)
		((ManagedMemory&)C[q]).read(fp);
	fclose(fp);
	#undef setFname
}

DataRptrCollection InvertChi::hessian(const DataRptrCollection& dV)
{	//Compute -dn = -chi * dV (first order perturbation theory)
	DataRptrCollection dn(dV.size());
	for(int q=0; q<e.eInfo.nStates; q++)
	{	int s = e.eInfo.qnums[q].index();
		for(unsigned b1=0; b1<eigs[q].size()-1; b1++)
		{	complexDataRptr conjIpsi1 = conj(I(C[q].getColumn(b1)));
			for(unsigned b2=b1+1; b2<eigs[q].size(); b2++)
				if(fabs(F[q][b1]-F[q][b2])>e.cntrl.occupiedThreshold)
				{	complexDataRptr n12 = conjIpsi1 * I(C[q].getColumn(b2));
					DataRptr n12real = Real(n12);
					DataRptr n12imag = Imag(n12);
					dn[s] += (e.eInfo.qnums[q].weight
							* (F[q][b1]-F[q][b2]) / (eigs[q][b2]-eigs[q][b1])
							* 2*e.gInfo.dV)
						* (n12real*dot(n12real,dV[s]) + n12imag*dot(n12imag,dV[s]));
				}
		}
	}
	return dn;
}

DataRptrCollection InvertChi::precondition(const DataRptrCollection& v)
{	return v;
}
