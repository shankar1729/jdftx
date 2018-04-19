/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#include <fluid/FluidSolverParams.h>
#include <electronic/SpeciesInfo.h>
#include <core/Units.h>

FluidSolverParams::FluidSolverParams()
: T(298*Kelvin), P(1.01325*Bar), epsBulkOverride(0.), epsInfOverride(0.), verboseLog(false), solveFrequency(FluidFreqDefault),
components(components_), solvents(solvents_), cations(cations_), anions(anions_),
vdwScale(0.75), pCavity(0.), lMax(3), cavityScale(1.), ionSpacing(0.),
zMask0(0.), zMaskH(0.), zMaskIonH(0.), zMaskSigma(0.5),
linearDielectric(false), linearScreening(false), nonlinearSCF(false), screenOverride(0.)
{
}

void FluidSolverParams::addComponent(const std::shared_ptr<FluidComponent>& component)
{	components_.push_back(component);
	switch(component->type)
	{	case FluidComponent::Solvent: solvents_.push_back(component); break;
		case FluidComponent::Cation: cations_.push_back(component); break;
		case FluidComponent::Anion: anions_.push_back(component); break;
	}
}


void FluidSolverParams::setCDFTparams()
{
	if(solvents.size() == 1)
	{
		switch(solvents[0]->name)
		{
			case FluidComponent::H2O:
				vdwScale = 0.540;
				break;
			case FluidComponent::CHCl3:
				vdwScale = 0.393;
				break;
			case FluidComponent::CCl4:
				vdwScale = 0.407;
				break;
			case FluidComponent::CH3CN:
			default:
			vdwScale = 0.488;
			initWarnings += "WARNING: Classical DFT has not been parameterized for this solvent, using 0.488 as the Van der Waals scale factor!\n";
			break;
				
		}
	}
	else
	{
	        vdwScale = 0.488;
		initWarnings += "WARNING: Classical DFT has not been parameterized for solvent mixtures, using 0.488 as the Van der Waals scale factor!\n";
	}
}

void FluidSolverParams::setPCMparams()
{
	assert(solvents.size()==1);
	const double dyn_per_cm = (1e-5*Newton)/(1e-2*meter);
	const double GPa = 1e9*Pascal;
	
	//Set common default values:
	cavityTension = 0.;
	cavityPressure = 0.;
	cavityScale = 1.;
	vdwScale = 1.;
	
	//Set PCM fit parameters:
	switch(pcmVariant)
	{	case PCM_SaLSA:
		{	nc = 1.42e-3;
			sigma = sqrt(0.5);
			switch(solvents[0]->name)
			{
				case FluidComponent::H2O:
					vdwScale = 0.50;
					break;
				case FluidComponent::CHCl3:
					vdwScale = 0.88;
					break;
				case FluidComponent::CCl4:
					vdwScale = 1.06;
					break;
				case FluidComponent::CH3CN:
					vdwScale = 0.37; //Tentative: using toluene solvation energy alone
					break;
				default:
					vdwScale = 1.;
					initWarnings += "WARNING: SALSA has not been parametrized for this solvent, using 1.0 as the Van der Waals scale factor!\n";
					break;
			}
			assert(fluidType == FluidSaLSA);
			break;
		}
		case PCM_CANDLE:
		{	nc = 1.42e-3;
			sigma = sqrt(0.5);
			switch(solvents[0]->name)
			{	case FluidComponent::CH3CN:
					Ztot = 16;
					eta_wDiel = 3.15;
					sqrtC6eff = 2.21;
					pCavity = -31.;
					break;
				case FluidComponent::Octanol:
					Ztot = 56;
					eta_wDiel = 5.5;
					sqrtC6eff = 16.;
					pCavity = 0;
					initWarnings += "WARNING: CANDLE parameters for Octanol have not yet been fit (initial guess only)\n";
					break;
				case FluidComponent::DMSO:
					Ztot = 26;
					eta_wDiel = 3.8;
					sqrtC6eff = 8.;
					pCavity = 20.;
					initWarnings += "WARNING: CANDLE parameters for DMSO have not yet been fit (initial guess only)\n";
					break;
				case FluidComponent::H2O:
				default:
					Ztot = 8;
					eta_wDiel = 1.46;
					sqrtC6eff = 0.770;
					pCavity = 36.5;
					if(solvents[0]->name != FluidComponent::H2O)
						initWarnings += "WARNING: CANDLE LinearPCM has not been parametrized for this solvent, using fit parameters for water\n";
					break;
			}
			if(fluidType != FluidLinearPCM) initWarnings += "WARNING: CANDLE has only been parametrized for LinearPCM.\n";
			break;
		}
		case PCM_SGA13:
		{	nc = 1e-2;
			sigma = 0.6;
			switch(solvents[0]->name)
			{
				case FluidComponent::H2O:
					vdwScale = 0.538;
					break;
				case FluidComponent::CHCl3:
					vdwScale = 1.315;
					break;
				case FluidComponent::CCl4:
					vdwScale = 1.238;
					break;
				case FluidComponent::CH3CN:
					vdwScale = 0.6; //Tentative: using toluene and self solvation energies alone
					break;
				default:
					vdwScale = 1.;
					initWarnings += "WARNING: SGA13 has not been parametrized for this solvent, using 1.0 as the Van der Waals scale factor!\n";
					break;
			}
			break;
		}
		case PCM_GLSSA13:
		{	
			switch(solvents[0]->name)
			{
				case FluidComponent::CHCl3:
				{	nc = 2.4e-05;
					sigma = 0.6;
					cavityTension = -9.066e-6;
					break;
				}
				case FluidComponent::CCl4:
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 1.15e-4;
							sigma = 0.6;
							cavityTension = -8.99e-06;
							break;
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with CCl4 as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
					break;
				}
				case FluidComponent::CH3CN:
				{	nc = 1.8e-4;
					sigma = 0.6;
					cavityTension = -6.3e-7;
					break;
				}
				case FluidComponent::CH2Cl2:
				{	nc = 9.3e-4;
					sigma = 0.6;
					cavityTension = -2.7e-6;
					break;
				}
				case FluidComponent::Ethanol:
				{	nc = 1.3e-3;
					sigma = 0.6;
					cavityTension = -5.1e-6;
					break;
				}
				case FluidComponent::Methanol:
				{	nc = 6.5e-4;
					sigma = 0.6;
					cavityTension = -5.2e-6;
					break;
				}
				case FluidComponent::EthylEther:
				{	nc = 2.63e-4;
					sigma = 0.6;
					cavityTension = -1.08e-5;
					break;
				}
				case FluidComponent::Chlorobenzene:
				{	nc = 4.28e-5;
					sigma = 0.6;
					cavityTension = -6.2e-06;
					break;
				}
				case FluidComponent::Isobutanol:
				{	nc = 1.51e-3;
					sigma = 0.6;
					cavityTension = -8.96e-06;
					break;
				}
				case FluidComponent::EC:
				{	nc = 1.8e-3;
					sigma = 0.6;
					cavityTension = 1.55e-5;
					break;
				}
				case FluidComponent::PC:
				{	nc = 9.8e-4;
					sigma = 0.6;
					cavityTension = 9.53e-6;
					break;
				}
				case FluidComponent::EthyleneGlycol:
				{	switch(fluidType)
					{	case FluidLinearPCM:
						{	nc = 5.4e-4;
							sigma = 0.6;
							cavityTension = 1.15e-5;
							break;
						}
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with Ethylene Glycol as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
				}
				case FluidComponent::Glyme:
				{	switch(fluidType)
					{	case FluidLinearPCM:
						{	nc = 8.36e-5;
							sigma = 0.6;
							cavityTension = -8.03e-06;
							break;
						}
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with Glyme as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
					break;
				}
				case FluidComponent::CarbonDisulfide:
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 1.6e-5;
							sigma = 0.6;
							cavityTension = -6.8e-06;
							break;
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with CarbonDisulfide as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
					break;
				}
				case FluidComponent::THF:
				{	nc = 1.6e-3;
					sigma = 0.6;
					cavityTension = -1.7e-06;
					break;
				}
				case FluidComponent::DMSO:
				{	nc = 9.5e-4;
					sigma = 0.6;
					cavityTension = 8.42e-06;
					break;
				}
				default: // For water and unparametrized fluids
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 3.7e-4;
							sigma = 0.6;
							cavityTension = 5.4e-6;
							break;
						case FluidNonlinearPCM:
							nc = 1.0e-3;
							sigma = 0.6;
							cavityTension = 9.5e-6;
							break;
						default: //Other fluids do not use these parameters
							break;
					}
					if(solvents[0]->name != FluidComponent::H2O)
					{	initWarnings +=
							"WARNING: PCM variant GLSSA has not been parametrized for this solvent; using bulk\n"
							"   surface tension as effective cavity tension and water parameters for cavity size.\n";
						cavityTension = solvents[0]->sigmaBulk;
					}
					break;
				}
			}
			break;
		}
		case PCM_LA12:
		{	nc = 7e-4;
			sigma = 0.6;
			if(fluidType == FluidNonlinearPCM)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has not been parametrized for NonlinearPCM; using LinearPCM fit parameters.\n";
			if( (fluidType==FluidLinearPCM || fluidType==FluidNonlinearPCM) && solvents[0]->name != FluidComponent::H2O)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has been fit only for H2O; using nc and sigma from H2O fit.\n";
			break;
		}
		case PCM_SoftSphere:
		{	sigma = 0.5; //in bohrs for this model
			switch(solvents[0]->name)
			{	case FluidComponent::Ethanol:
				{	cavityTension = -4.0*dyn_per_cm;
					cavityScale = 1.22;
					break;
				}
				default: // For water and unparametrized fluids
				{	if(fluidType==FluidLinearPCM)
					{	cavityTension = 9.67e-6;
						cavityScale = 1.11;
					}
					else
					{	cavityTension = 1.02e-5;
						cavityScale = 1.00;
					}
					//Warn when not explicitly parametrized:
					if(solvents[0]->name != FluidComponent::H2O)
					{	initWarnings +=
							"WARNING: PCM variant SoftSphere has not been parametrized for this solvent; using bulk\n"
							"   surface tension as effective cavity tension and water parameters for everything else.\n";
						cavityTension = solvents[0]->sigmaBulk;
					}
				}
			}
		}
		case_PCM_SCCS_any:
		{	if(fluidType != FluidLinearPCM) initWarnings += "WARNING: SCCS has only been parametrized for LinearPCM.\n";
			if(solvents[0]->name != FluidComponent::H2O)
			{	initWarnings += 
					"WARNING: SCCS variants have not been parametrized for this solvent; using water parameters\n";
			}
			rhoDelta = 1e-4; //common to all variants
			switch(pcmVariant)
			{	case PCM_SCCS_g09:      rhoMin=1.00e-4; rhoMax=1.50e-3; cavityTension=2.50*dyn_per_cm; break;
				case PCM_SCCS_g03:      rhoMin=1.00e-4; rhoMax=5.00e-3; cavityTension=11.5*dyn_per_cm; break;
				case PCM_SCCS_g03p:     rhoMin=3.00e-4; rhoMax=3.00e-3; cavityTension=12.0*dyn_per_cm; break;
				case PCM_SCCS_g09beta:  rhoMin=1.00e-4; rhoMax=1.50e-3; cavityTension=11.0*dyn_per_cm; cavityPressure=-0.08*GPa; break;
				case PCM_SCCS_g03beta:  rhoMin=1.00e-4; rhoMax=5.00e-3; cavityTension=50.0*dyn_per_cm; cavityPressure=-0.35*GPa; break;
				case PCM_SCCS_g03pbeta: rhoMin=3.00e-4; rhoMax=3.00e-3; cavityTension=20.0*dyn_per_cm; cavityPressure=-0.08*GPa; break;
				case PCM_SCCS_cation:   rhoMin=2.00e-4; rhoMax=3.50e-3; cavityTension=5.00*dyn_per_cm; cavityPressure=+0.125*GPa; break;
				case PCM_SCCS_anion:    rhoMin=2.40e-3; rhoMax=1.55e-2; cavityTension=0.00*dyn_per_cm; cavityPressure=+0.450*GPa; break;
				default: break;
			}
			break;
		}
	}
}

double FluidSolverParams::getAtomicRadius(const SpeciesInfo& sp) const
{	static const int Zmax = 103;
	static double uffDiaAngst[Zmax] = {
		2.886, 2.362, 2.451, 2.745, 4.083, 3.851, 3.100, 3.500, 3.364, 3.243, //Z -> 10 (Z=7 i.e. N overridden to Bondi vdW radius)
		2.983, 3.021, 4.499, 4.295, 4.147, 4.035, 3.947, 3.868, 3.812, 3.399, //Z -> 20
		3.295, 3.175, 3.144, 3.023, 2.961, 2.912, 2.872, 2.834, 3.495, 2.763, //Z -> 30
		4.383, 4.280, 4.230, 4.205, 4.189, 4.141, 4.114, 3.641, 3.345, 3.124, //Z -> 40
		3.165, 3.052, 2.998, 2.963, 2.929, 2.899, 3.148, 2.848, 4.463, 4.392, //Z -> 50
		4.420, 4.470, 4.500, 4.404, 4.517, 3.703, 3.522, 3.556, 3.606, 3.575, //Z -> 60
		3.547, 3.520, 3.493, 3.368, 3.451, 3.428, 3.409, 3.391, 3.374, 3.355, //Z -> 70
		3.640, 3.141, 3.170, 3.069, 2.954, 3.120, 2.840, 2.754, 3.293, 2.705, //Z -> 80
		4.347, 4.297, 4.370, 4.709, 4.750, 4.765, 4.900, 3.677, 3.478, 3.396, //Z -> 90
		3.424, 3.395, 3.424, 3.424, 3.381, 3.326, 3.339, 3.313, 3.299, 3.286, //Z -> 100
		3.274, 3.248, 3.236 }; //Z -> 103
	if(sp.atomicRadiusOverride) return *(sp.atomicRadiusOverride); //Return override value, if available
	assert(sp.atomicNumber && sp.atomicNumber <= Zmax);
	return uffDiaAngst[sp.atomicNumber-1]*0.5*Angstrom; //convert to radius in bohrs
}


bool FluidSolverParams::needsVDW() const
{	switch(fluidType)
	{	case FluidNone:
			return false;
		case FluidLinearPCM:
		case FluidNonlinearPCM:
			return (pcmVariant==PCM_SGA13 || pcmVariant==PCM_CANDLE);
		case FluidSaLSA:
		case FluidClassicalDFT:
		default:
			return true;
	}
}

bool FluidSolverParams::ionicScreening() const
{	return cations.size() && anions.size();
}
