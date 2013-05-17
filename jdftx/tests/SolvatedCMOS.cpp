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
#include <fluid/IdealGasMuEps.h>
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/IdealGasPomega.h>
#include <fluid/IdealGasMonoatomic.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_LJ.h>

const double epsSi = 11.8;
const double epsSiO2 = 3.9;

const double cm = 1e-2*meter;
const double cm3 = pow(cm,3);

//Semiconductor properties (denisties are per electrical distance, hence scaled by dielectric constant):
const double Ndope = 3e17/cm3 * epsSi; //dopant density
static enum { Ptype, Ntype } dopeType = Ptype; //doping type
const double Nintrinsic = 1e14/cm3 * epsSi; //equlibirum electron-hole densities in intrinsic silicon
const double Vbuiltin = 0;//0.8*eV; //built in potential between silicon and floating gate (positive puts holes on floating gate)

const double sigmaFG = 3e11/pow(cm,2); //excess electrons per unit area on floating gate

//Electrical distances (actual length/dielectric constant)
const double Lsilicon = 1000*Angstrom / epsSi;
const double Lox1 = 100*Angstrom / epsSiO2; //between silicon and floating gate
const double Lox2 = 350*Angstrom / epsSiO2; //between floating gate and liquid interface

//Actual distance (dielectric response handled microscopically)
const double Lliq = 100*Angstrom;

//Compute the various interface distances from cell center based on the above:
const double z_Si_Ox1 = Lsilicon;
const double z_Ox1_Ox2 = z_Si_Ox1 + Lox1;
const double z_Ox2_liq = z_Ox1_Ox2 + Lox2;
const double zEnd = z_Ox2_liq + Lliq;

const int nPoints = 2048;
const double hGrid = (2*zEnd)/nPoints;

//Set the repulsive potentials (of strength Vrep) that keep the carriers in Si and the ions in the liquid etc.
void setWalls(int i, vector3<> r, double Vrep,
	double* VO, double* VH, double* Vcation, double* Vanion, double* Ve, double* Vh)
{	double z = fabs(r[2] - zEnd); //symmetric about center of cell (to prevent net dipole)
	VO[i]=VH[i]=Vcation[i]=Vanion[i] = z>z_Ox2_liq ? 0 : Vrep; //water and ions only in liquid area
	Ve[i]=Vh[i] = z<z_Si_Ox1 ? 0 : Vrep; //electrons and holes only in silicon
}

//Set the background charge density in semiconductor due to dopant ions
void setDopantBG(int i, vector3<> r, double rho, double* rhoDopantBG)
{	double z = fabs(r[2] - zEnd); //symmetric about center of cell (to prevent net dipole)
	rhoDopantBG[i] = z<z_Si_Ox1 ? rho : 0;
}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	printf("Grid spacing = %lf A\n", hGrid/Angstrom);
	GridInfo gInfo;
	gInfo.S = vector3<int>(1, 1, nPoints);
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();

	double T = 298*Kelvin;
	FluidComponent componentH2O(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
	FluidComponent componentCation(FluidComponent::Sodium, T, FluidComponent::MeanFieldLJ);
	FluidComponent componentAnion(FluidComponent::Chloride, T, FluidComponent::MeanFieldLJ);
	componentCation.Nbulk = 0.02*mol/liter;
	componentAnion.Nbulk = 0.02*mol/liter;
	FluidComponent componentHoles(FluidComponent::CustomCation, T, FluidComponent::MeanFieldLJ);
	FluidComponent componentElectrons(FluidComponent::CustomAnion, T, FluidComponent::MeanFieldLJ);
	componentHoles.molecule.setModelMonoatomic("h+", -1., 0*Angstrom);
	componentElectrons.molecule.setModelMonoatomic("e-", +1., 0*Angstrom);

	FluidMixture fluidMixture(gInfo, T);
	componentH2O.addToFluidMixture(&fluidMixture);
	componentCation.addToFluidMixture(&fluidMixture);
	componentAnion.addToFluidMixture(&fluidMixture);
	componentHoles.addToFluidMixture(&fluidMixture);
	componentElectrons.addToFluidMixture(&fluidMixture);
	//Mixing functionals:
	Fmix_LJ fsolvCation(&fluidMixture, &componentH2O, &componentCation, 2e-3, 3.0*Angstrom);
	Fmix_LJ fsolvAnion(&fluidMixture, &componentH2O, &componentAnion, 2e-3, 3.0*Angstrom);
	
	double p = 1.01325*Bar;
	printf("pV = %le\n", p*gInfo.detR);
	fluidMixture.initialize(p);

	//Set semiconductor properties:
	DataRptr rhoDopantBG, rhoBuiltin;
	{	printf("Initializing semiconductor:\n");
		//Figure out bulk concentrations of majority and minority carriers:
		double Nmajor = 0.5*(fabs(Ndope) +  sqrt(pow(Ndope,2) + 4*pow(Nintrinsic,2)));
		double Nminor = pow(Nintrinsic,2)/Nmajor;
		double Nelectrons, Nholes, rhoBG;
		switch(dopeType)
		{	case Ptype: Nelectrons=Nminor; Nholes=Nmajor; rhoBG=+fabs(Ndope); break;
			case Ntype: Nelectrons=Nmajor; Nholes=Nminor; rhoBG=-fabs(Ndope); break;
		}
		printf("\tBulk electron density = %le/cm3\n", Nelectrons*cm3);
		printf("\tBulk hole density = %le/cm3\n", Nholes*cm3);
		printf("\tDopant ionic charge density = %le/cm3\n", rhoBG*cm3);
		componentHoles.idealGas->overrideBulk(Nholes, 0);
		componentElectrons.idealGas->overrideBulk(Nelectrons, 0);
		nullToZero(rhoDopantBG, gInfo);
		applyFunc_r(gInfo, setDopantBG, rhoBG, rhoDopantBG->data());

		//Handle built-in potential:
		double sigmaBuiltin = Vbuiltin/(4*M_PI*Lox1);
		DataRptr rhoSheet; initZero(rhoSheet, gInfo);
		rhoSheet->data()[0] = sigmaBuiltin/hGrid;
		RealKernel gauss(gInfo); initGaussianKernel(gauss, 0.5*Angstrom);
		DataGptr transFG1(DataG::alloc(gInfo)); initTranslation(transFG1, vector3<>(0,0,zEnd+z_Ox1_Ox2));
		DataGptr transFG2(DataG::alloc(gInfo)); initTranslation(transFG2, vector3<>(0,0,zEnd-z_Ox1_Ox2));
		DataGptr transSubs1(DataG::alloc(gInfo)); initTranslation(transSubs1, vector3<>(0,0,zEnd+z_Si_Ox1));
		DataGptr transSubs2(DataG::alloc(gInfo)); initTranslation(transSubs2, vector3<>(0,0,zEnd-z_Si_Ox1));
		rhoBuiltin = I((transSubs1*0.0+transSubs2*0.0-transFG1-transFG2)*(gauss*J(rhoSheet)));
	}

	//----- Repulsive potentials to keep things in place -----
	nullToZero(componentH2O.idealGas->V, gInfo);
	nullToZero(componentCation.idealGas->V, gInfo);
	nullToZero(componentAnion.idealGas->V, gInfo);
	nullToZero(componentHoles.idealGas->V, gInfo);
	nullToZero(componentElectrons.idealGas->V, gInfo);
	applyFunc_r(gInfo, setWalls, 1.0,
		componentH2O.idealGas->V[0]->data(), componentH2O.idealGas->V[1]->data(),
		componentCation.idealGas->V[0]->data(), componentAnion.idealGas->V[0]->data(),
		componentHoles.idealGas->V[0]->data(), componentElectrons.idealGas->V[0]->data() );

	//----- Initialize external charge -----
	DataRptr rhoFloatingGate;
	{	DataRptr rhoFGsingle; initZero(rhoFGsingle, gInfo);
		rhoFGsingle->data()[0] = sigmaFG/hGrid; //sheet charge at z=0
		RealKernel gauss(gInfo); initGaussianKernel(gauss, 0.5*Angstrom);
		DataGptr trans1(DataG::alloc(gInfo)); initTranslation(trans1, vector3<>(0,0,zEnd+z_Ox1_Ox2));
		DataGptr trans2(DataG::alloc(gInfo)); initTranslation(trans2, vector3<>(0,0,zEnd-z_Ox1_Ox2));
		rhoFloatingGate = I(gauss*(trans1+trans2)*J(rhoFGsingle));
	}

	//Put together all the external charges:
	fluidMixture.rhoExternal = J(rhoFloatingGate + rhoDopantBG + rhoBuiltin);

	//---- G=0 constraint -----
	//fluidMixture.d0calc = zeroCenter;

	//----- Initialize state -----
	fluidMixture.initState(0.15);

	//----- FDtest and CG -----
	MinimizeParams mp;
	mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
	mp.nIterations=1000;
	mp.knormThreshold=1e-14*pow(hGrid,2);
	mp.fdTest = true;
	
	puts("Starting CG:");
	TIME("minimize", stdout,
		fluidMixture.minimize(mp);
	);

	DataRptrCollection N; DataGptr grad_rhoExternalTilde;
	TIME("getFreeEnergy", stdout,
		fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, 0, &grad_rhoExternalTilde));
	);

	//Scale densities by bulk values:
	const std::vector<const FluidComponent*>& component = fluidMixture.getComponents();
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *(component[ic]);
		double Nbulk = c.idealGas->get_Nbulk();
		for(unsigned j=0; j<c.molecule.sites.size(); j++)
			N[c.offsetDensity+j] *= 1./(Nbulk * c.molecule.sites[j]->positions.size());
	}

	DataRptr dtot = I(grad_rhoExternalTilde - 4*M_PI*Linv(O(fluidMixture.rhoExternal)));
	DataRptr dtotNoGzero = dtot - integral(dtot)/gInfo.detR;
	FILE* fp = fopen("SolvatedCMOS.N", "w");
	for(int i=0; i<gInfo.nr; i++)
	{	fprintf(fp, "%le", hGrid*i/Angstrom);
		for(unsigned j=0; j<N.size(); j++)
			fprintf(fp, "\t%le", N[j]->data()[i]);
		fprintf(fp, "\t%le\t%le", dtot->data()[i]/eV, dtotNoGzero->data()[i]/eV);
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	int i_Si_Ox1 = ceil((zEnd + z_Si_Ox1)/hGrid)+2;
	int i_Ox1_Ox2 = floor((zEnd+z_Ox1_Ox2)/hGrid)-2;
	double Eox1 = (dtot->data()[i_Ox1_Ox2] - dtot->data()[i_Si_Ox1])/((eV/Angstrom)*hGrid*(i_Ox1_Ox2-i_Si_Ox1));
	printf("Eox1 = %le V/A\n", Eox1);
}
