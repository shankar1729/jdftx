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

#include <fluid/FluidMixture.h>
#include <core/ScalarFieldIO.h>
#include <core/string.h>
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>

void initCavity(size_t i,double nc, const double* nEl, double* phiO, double* phiH)
{	if(nEl[i]>nc) { phiO[i]=1.0; phiH[i]=0.0; }
	else  { phiO[i]=0.0; phiH[i]=0.0; }
}

int main(int argc, char** argv)
{	initSystem(argc, argv);
	if(argc != 3)
	{	logPrintf("\tUsage: CavitationH2O <output-file> <binary-density-file>\n"
			"\t\t<output-file>: JDFTx output file (used to get lattice vectors and sample counts)\n"
			"\t\t<binary-density-file>: Electron-density (JDFTx format) which is thresholded to form cavity\n");
		return 1;
	}
	initSystem(argc, argv);
	
	ifstream ifs(argv[1]);
	if(!ifs.is_open()) { logPrintf("Could not open '%s' for reading\n", argv[1]); return 1; }

	GridInfo gInfo; double pcm_nc, pcm_sigma;

	bool gotR=false, gotS=false, got_nc=true;
	while(!ifs.eof())
	{	string line; getline(ifs, line);
		if(line.substr(0,11)=="jdft1-shape") //TODO: change command name to pcm-params
		{	istringstream ss(line.substr(11));
			ss >> pcm_nc >> pcm_sigma;
			got_nc=true;
		}
		if(line=="R = ")
		{	string lBkt, rBkt;
			ifs >> lBkt >>  gInfo.R(0,0) >> gInfo.R(0,1) >> gInfo.R(0,2) >> rBkt;
			ifs >> lBkt >>  gInfo.R(1,0) >> gInfo.R(1,1) >> gInfo.R(1,2) >> rBkt;
			ifs >> lBkt >>  gInfo.R(2,0) >> gInfo.R(2,1) >> gInfo.R(2,2) >> rBkt;
			gotR=true;
		}
		if(line.substr(0,25)=="Chosen fftbox size, S = [")
		{	istringstream ss(line.substr(25));
			ss >> gInfo.S[0] >> gInfo.S[1] >> gInfo.S[2];
			gotS=true;
		}
	}
	if(!got_nc) logPrintf("Could not read critical density nc for PCM shape from '%s'\n", argv[1]);
	if(!gotR) logPrintf("Could not read lattice vectors (R[][]) from '%s'\n", argv[1]);
	if(!gotS) logPrintf("Could not read sample counts (S[]) from '%s'\n", argv[1]);
	if(!(got_nc && gotR && gotS)) return 1;

	logPrintf("\nParsed '%s' and found (vectors in rows):\n", argv[1]);
	logPrintf("\tR[0] = [ %lf\t%lf\t%lf ];\n",gInfo.R(0,0), gInfo.R(1,0), gInfo.R(2,0));
	logPrintf("\tR[1] = [ %lf\t%lf\t%lf ];\n",gInfo.R(0,1), gInfo.R(1,1), gInfo.R(2,1));
	logPrintf("\tR[2] = [ %lf\t%lf\t%lf ];\n",gInfo.R(0,2), gInfo.R(1,2), gInfo.R(2,2));
	logPrintf("\tS = [ %d %d %d ];\n\n", gInfo.S[0], gInfo.S[1], gInfo.S[2]);
	gInfo.initialize();


	double T = 298*Kelvin;
	FluidComponent component(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
	component.s2quadType = Quad7design_24;
	
	FluidMixture fluidMixture(gInfo, T);
	component.addToFluidMixture(&fluidMixture);
	double p = 1.01325*Bar;
	logPrintf("pV = %le\n", p*gInfo.detR);
	fluidMixture.initialize(p);

	//----- Initialize external potential -----
	logPrintf("\nReading electron density from '%s' ... ", argv[1]); fflush(stdout);
	ScalarField nElectronic(ScalarFieldData::alloc(gInfo));
	loadRawBinary(nElectronic, argv[2]);
	logPrintf("SumCheck: %lf electrons\n\n", integral(nElectronic));

	logPrintf("Initializing cavity with nc = %le ... ", pcm_nc); fflush(stdout);
	IdealGas& idgas = *(component.idealGas);
	nullToZero(idgas.V, gInfo);
	threadedLoop(initCavity, gInfo.nr, pcm_nc, nElectronic->data(), idgas.V[0]->data(), idgas.V[1]->data());
	logPrintf("CavityVolume: %lf bohr^3\n\n", integral(idgas.V[0]));


	RealKernel gaussian(gInfo);
	initGaussianKernel(gaussian, sqrt(gInfo.h[0].length_squared()+gInfo.h[1].length_squared()+gInfo.h[2].length_squared()));

	for (int i=0; i<2; i++) idgas.V[i] = 0.005 * I(gaussian * J(idgas.V[i]));
	fluidMixture.initState(0.15);

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = 2*gInfo.nr;
	mp.nIterations=200;
	mp.knormThreshold=1e-11;
	mp.dirUpdateScheme=MinimizeParams::HestenesStiefel;
	
	logPrintf("Starting CG:\n");
	TIME("minimize", globalLog,
		fluidMixture.minimize(mp);
	);
	
	ScalarFieldArray N(2);
	double phiFinal = fluidMixture.getFreeEnergy(&N);
	logPrintf("\nCavitation energy = %le Eh\n\n", phiFinal);
	
	finalizeSystem();
	return 0;
}
