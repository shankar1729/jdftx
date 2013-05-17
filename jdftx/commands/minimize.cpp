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

#include <commands/command.h>
#include <electronic/Everything.h>

EnumStringMap<MinimizeParams::DirectionUpdateScheme> dirUpdateMap
(	MinimizeParams::PolakRibiere, "PolakRibiere",
	MinimizeParams::FletcherReeves, "FletcherReeves",
	MinimizeParams::HestenesStiefel, "HestenesStiefel",
	MinimizeParams::SteepestDescent, "SteepestDescent"
);

EnumStringMap<MinimizeParams::LinminMethod> linminMap
(	MinimizeParams::Relax, "Relax",
	MinimizeParams::Quad, "Quad",
	MinimizeParams::Cubic, "Cubic"
);

//An enum entry for each configurable option of MinimizeParams
enum MinimizeParamsMember
{	MPM_dirUpdateScheme,
	MPM_linminMethod,
	MPM_nIterations,
	MPM_knormThreshold,
	MPM_energyDiffThreshold,
	MPM_nEnergyDiff,
	MPM_nDirResetNum,
	MPM_nDirResetDen,
	MPM_alphaTstart,
	MPM_alphaTmin,
	MPM_updateTestStepSize,
	MPM_alphaTreduceFactor,
	MPM_alphaTincreaseFactor,
	MPM_nAlphaAdjustMax,
	MPM_fdTest,
	MPM_Delim //!< delimiter to detect end of input
};

EnumStringMap<MinimizeParamsMember> mpmMap
(	MPM_dirUpdateScheme, "dirUpdateScheme",
	MPM_linminMethod, "linminMethod",
	MPM_nIterations, "nIterations",
	MPM_knormThreshold, "knormThreshold",
	MPM_energyDiffThreshold, "energyDiffThreshold",
	MPM_nEnergyDiff, "nEnergyDiff",
	MPM_nDirResetNum, "nDirResetNum",
	MPM_nDirResetDen, "nDirResetDen",
	MPM_alphaTstart, "alphaTstart",
	MPM_alphaTmin, "alphaTmin",
	MPM_updateTestStepSize, "updateTestStepSize",
	MPM_alphaTreduceFactor, "alphaTreduceFactor",
	MPM_alphaTincreaseFactor, "alphaTincreaseFactor",
	MPM_fdTest, "fdTest",
	MPM_nAlphaAdjustMax, "nAlphaAdjustMax"
);
EnumStringMap<MinimizeParamsMember> mpmDescMap
(	MPM_dirUpdateScheme, dirUpdateMap.optionList() + " (search direction method)",
	MPM_linminMethod, linminMap.optionList() + " (line minimization method)",
	MPM_nIterations, "maximum iterations (single point calculation if 0)",
	MPM_knormThreshold, "convergence threshold for gradient (preconditioned) norm",
	MPM_energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
	MPM_nEnergyDiff, "number of iteration pairs that must satisfy energyDiffThreshold",
	MPM_nDirResetDen, "<int>  / reset search direction when linmin fails <nDirResetNum> \\",
	MPM_nDirResetNum, "<int>  \\     out of <nDirResetDen> successive iterations         /",
	MPM_alphaTstart, "initial test step size (constant step-size factor for Relax linmin)",
	MPM_alphaTmin, "minimum test step size",
	MPM_updateTestStepSize, boolMap.optionList() + ", whether test step size is updated",
	MPM_alphaTreduceFactor, "step size reduction factor when energy increases in linmin",
	MPM_alphaTincreaseFactor, "step size increase factor when curvature is wrong in linmin",
	MPM_nAlphaAdjustMax, "maximum step-size adjustments per linmin",
	MPM_fdTest, boolMap.optionList() + ", whether to perform a finite difference test"
);


//! Minimization parameters for ions, electrons or fluid
struct CommandMinimize : public Command
{
	CommandMinimize(string systemName) : Command(systemName + "-minimize")
	{
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "where possible keys and value types are:"
			+ addDescriptions(mpmMap.optionList(), linkDescription(mpmMap, mpmDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	MinimizeParams& mp = target(e);
		while(true)
		{	MinimizeParamsMember key;
			pl.get(key, MPM_Delim, mpmMap, "key");
			switch(key)
			{	case MPM_dirUpdateScheme: pl.get(mp.dirUpdateScheme, MinimizeParams::PolakRibiere, dirUpdateMap, "dirUpdateScheme", true); break;
				case MPM_linminMethod: pl.get(mp.linminMethod, MinimizeParams::Quad, linminMap, "linminMethod", true); break;
				case MPM_nIterations: pl.get(mp.nIterations, 0, "nIterations", true); break;
				case MPM_knormThreshold: pl.get(mp.knormThreshold, 0., "knormThreshold", true); break;
				case MPM_energyDiffThreshold: pl.get(mp.energyDiffThreshold, 0., "energyDiffThreshold", true); break;
				case MPM_nEnergyDiff: pl.get(mp.nEnergyDiff, 0, "nEnergyDiff", true); break;
				case MPM_nDirResetNum: pl.get(mp.nDirResetNum, 0, "nDirResetNum", true); break;
				case MPM_nDirResetDen: pl.get(mp.nDirResetDen, 0, "nDirResetDen", true); break;
				case MPM_alphaTstart: pl.get(mp.alphaTstart, 0., "alphaTstart", true); break;
				case MPM_alphaTmin: pl.get(mp.alphaTmin, 0., "alphaTmin", true); break;
				case MPM_updateTestStepSize: pl.get(mp.updateTestStepSize, true, boolMap, "updateTestStepSize", true); break;
				case MPM_alphaTreduceFactor: pl.get(mp.alphaTreduceFactor, 0., "alphaTreduceFactor", true); break;
				case MPM_alphaTincreaseFactor: pl.get(mp.alphaTincreaseFactor, 0., "alphaTincreaseFactor", true); break;
				case MPM_nAlphaAdjustMax: pl.get(mp.nAlphaAdjustMax, 0, "nAlphaAdjustMax", true); break;
				case MPM_fdTest: pl.get(mp.fdTest, false, boolMap, "fdTest", true); break;
				case MPM_Delim: return; //end of input
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	MinimizeParams& mp = target(e);
		logPrintf(" \\\n\tdirUpdateScheme      %s", dirUpdateMap.getString(mp.dirUpdateScheme));
		logPrintf(" \\\n\tlinminMethod         %s", linminMap.getString(mp.linminMethod));
		logPrintf(" \\\n\tnIterations          %d", mp.nIterations);
		logPrintf(" \\\n\tknormThreshold       %lg", mp.knormThreshold);
		logPrintf(" \\\n\tenergyDiffThreshold  %lg", mp.energyDiffThreshold);
		logPrintf(" \\\n\tnEnergyDiff          %d", mp.nEnergyDiff);
		logPrintf(" \\\n\tnDirResetNum         %d", mp.nDirResetNum);
		logPrintf(" \\\n\tnDirResetDen         %d", mp.nDirResetDen);
		logPrintf(" \\\n\talphaTstart          %lg", mp.alphaTstart);
		logPrintf(" \\\n\talphaTmin            %lg", mp.alphaTmin);
		logPrintf(" \\\n\tupdateTestStepSize   %s", boolMap.getString(mp.updateTestStepSize));
		logPrintf(" \\\n\talphaTreduceFactor   %lg", mp.alphaTreduceFactor);
		logPrintf(" \\\n\talphaTincreaseFactor %lg", mp.alphaTincreaseFactor);
		logPrintf(" \\\n\tnAlphaAdjustMax      %d", mp.nAlphaAdjustMax);
		logPrintf(" \\\n\tfdTest               %s", boolMap.getString(mp.fdTest));
	}

protected:
	virtual MinimizeParams& target(Everything&)=0; //!< derived class specifies where the parameters are stored
};


struct CommandElectronicMinimize : public CommandMinimize
{	CommandElectronicMinimize() : CommandMinimize("electronic") {}
    MinimizeParams& target(Everything& e) { return e.elecMinParams; }
    void process(ParamList& pl, Everything& e)
	{	//Use default value of 100 iterations from MinimizeParams.h
		e.elecMinParams.energyDiffThreshold = 1e-8;
		CommandMinimize::process(pl, e);
	}
}
commandElectronicMinimize;

struct CommandIonicMinimize : public CommandMinimize
{	CommandIonicMinimize() : CommandMinimize("ionic")
	{	emptyParamError =
			"   Note: nIterations defaults to 0 for ionic minimization,\n"
			"      and must be set manually to enable this feature.";
	}
    MinimizeParams& target(Everything& e) { return e.ionicMinParams; }
    void process(ParamList& pl, Everything& e)
	{	e.ionicMinParams.nIterations = 0; //override default value (100) in MinimizeParams.h
		e.ionicMinParams.energyDiffThreshold = 1e-6;
		e.ionicMinParams.knormThreshold = 1e-4;
		CommandMinimize::process(pl, e);
	}
}
commandIonicMinimize;

struct CommandFluidMinimize : public CommandMinimize
{	CommandFluidMinimize() : CommandMinimize("fluid")
	{	require("fluid");
	}
    MinimizeParams& target(Everything& e) { return e.fluidMinParams; }
    void process(ParamList& pl, Everything& e)
	{	const FluidType& fluidType = e.eVars.fluidParams.fluidType;
		switch(fluidType)
		{	case FluidLinearPCM:
			case FluidNonlocalPCM:
				e.fluidMinParams.nIterations = 400; //override default value (100) in MinimizeParams.h
				e.fluidMinParams.knormThreshold = (fluidType==FluidNonlocalPCM) ? 1e-12 : 1e-11;
				break;
			default:;
		}
		CommandMinimize::process(pl, e);
	}
}
commandFluidMinimize;

struct CommandLatticeMinimize : public CommandMinimize
{	CommandLatticeMinimize() : CommandMinimize("lattice")
	{	emptyParamError =
			"   Note: nIterations defaults to 0 for lattice minimization,\n"
			"      and must be set manually to enable this feature.";
	}
    MinimizeParams& target(Everything& e) { return e.latticeMinParams; }
    void process(ParamList& pl, Everything& e)
	{	e.latticeMinParams.nIterations = 0; //override default value (100) in MinimizeParams.h
		e.latticeMinParams.energyDiffThreshold = 1e-6;
		CommandMinimize::process(pl, e);
	}
}
commandLatticeMinimize;

struct CommandWannierMinimize : public CommandMinimize
{	CommandWannierMinimize() : CommandMinimize("wannier") {}
    MinimizeParams& target(Everything& e) { return e.dump.wannier.minParams; }
    void process(ParamList& pl, Everything& e)
	{	e.dump.wannier.minParams.energyDiffThreshold = 1e-8; //override default value (0.) in MinimizeParams.h
		CommandMinimize::process(pl, e);
	}
}
commandWannierMinimize;

struct CommandInverseKohnShamMinimize : public CommandMinimize
{	CommandInverseKohnShamMinimize() : CommandMinimize("inverseKohnSham") {}
    MinimizeParams& target(Everything& e) { return e.inverseKSminParams; }
    void process(ParamList& pl, Everything& e)
	{	e.inverseKSminParams.energyDiffThreshold = 1e-8; //override default value (0.) in MinimizeParams.h
		CommandMinimize::process(pl, e);
	}
}
commandInverseKohnShamMinimize;
