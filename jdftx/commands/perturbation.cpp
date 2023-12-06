/*-------------------------------------------------------------------
Copyright 2022 Brandon Li

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

#include <commands/minimize.h>
#include <electronic/PerturbationParams.h>

//An enum entry for each configurable option of MinimizeParams
enum VPTParamsMember
{	VPTp_nIterations,
	VPTp_algorithm,
	VPTp_residualTol,
	VPTp_residualDiffThreshold,
	VPTp_CGBypass,
	VPTp_recomputeResidual,
	VPTp_Delim //!< delimiter to detect end of input
};

EnumStringMap<VPTParamsMember> VPTpMap
(	VPTp_nIterations, "nIterations",
	VPTp_algorithm, "algorithm",
	VPTp_residualTol, "residualTol",
	VPTp_residualDiffThreshold, "residualDiffThreshold",
	VPTp_CGBypass, "CGBypass",
	VPTp_recomputeResidual, "recomputeResidual"
);
EnumStringMap<VPTParamsMember> VPTpDescMap
(	VPTp_nIterations, "maximum number of iterations",
	VPTp_algorithm, "CG|MINRES (linear solver algorithm)",
	VPTp_residualTol, "relative error of residual required to terminate solver (default value of 1e-4 is sufficient most of the time)",
	VPTp_residualDiffThreshold, "terminate if difference in consecutive residuals falls below threshold (MINRES only)",
	VPTp_CGBypass, "bypass CG stopping criterion (residual increases 3 iterations in a row)",
	VPTp_recomputeResidual, "recompute residual every ten iterations (CG only)"
);

EnumStringMap<PerturbationParams::Algorithms> algMap
(	PerturbationParams::CG, "CG",
	PerturbationParams::MINRES, "MINRES"
);

struct CommandPerturbMinimize : public Command
{
	CommandPerturbMinimize() : Command("perturb-minimize", "jdftx/Perturb")
	{	emptyParamError =
			"   Note: nIterations defaults to 0 and must be set manually to enable variational perturbation solver.\n";
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "where possible keys and value types are:"
			+ addDescriptions(VPTpMap.optionList(), linkDescription(VPTpMap, VPTpDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = true;
	}
	
    void process(ParamList& pl, Everything& e)
    {
		while(true)
		{	VPTParamsMember key;
			pl.get(key, VPTp_Delim, VPTpMap, "key");
			switch(key)
			{	case VPTp_nIterations: pl.get(e.vptParams.nIterations, 0, "nIterations", true); break;
				case VPTp_algorithm: pl.get(e.vptParams.algorithm, PerturbationParams::MINRES, algMap, "algorithm", true); break;
				case VPTp_residualTol: pl.get(e.vptParams.residualTol, 1e-4, "residualTol", true); break;
				case VPTp_residualDiffThreshold: pl.get(e.vptParams.residualDiffThreshold, 1e-4, "residualDiffThreshold", true); break;
				case VPTp_CGBypass: pl.get(e.vptParams.CGBypass, false, boolMap, "CGBypass", true); break;
				case VPTp_recomputeResidual: pl.get(e.vptParams.recomputeResidual, false, boolMap, "recomputeResidual", true); break;
				case VPTp_Delim: return; //end of input
			}
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	PerturbationParams& p = e.vptParams;
		logPrintf(" \\\n\tnIterations            %d", p.nIterations);
		logPrintf(" \\\n\talgorithm              %s", algMap.getString(p.algorithm));
		logPrintf(" \\\n\tresidualTol            %g", p.residualTol);
		logPrintf(" \\\n\tresidualDiffThreshold  %g", p.residualDiffThreshold);
		logPrintf(" \\\n\tCGBypass               %s", boolMap.getString(p.CGBypass));
		logPrintf(" \\\n\trecomputeResidual      %s", boolMap.getString(p.recomputeResidual));
	}
}
commandPerturbMinimize;


struct CommandPerturbVexternal : public Command
{
	CommandPerturbVexternal() : Command("perturb-Vexternal", "jdftx/Perturb")
	{
		format = "<filename>";
		comments =
			"Introduce external potential perturbation (in Hartrees).\n"
			"Must be specified as a real space field, or complex field for incommensurate perturbations.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.vptInfo.dVext = std::make_shared<VextPerturbation>(e);
		e.vptInfo.dVext->dVexternalFilename.resize(1);
		pl.get(e.vptInfo.dVext->dVexternalFilename[0], string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{	for(const string& filename: e.vptInfo.dVext->dVexternalFilename)
			logPrintf("%s ", filename.c_str());
	}
}
commandPerturbVexternal;


struct CommandPerturbRhoExternal : public Command
{
	CommandPerturbRhoExternal() : Command("perturb-rhoExternal", "jdftx/Perturb")
	{
		format = "<filename>";
		comments =
			"Introduce external charge density perturbation.\n"
			"Must be specified as a real space field, or complex field for incommensurate perturbations.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.vptInfo.drhoExt = std::make_shared<RhoPerturbation>(e);
		e.vptInfo.drhoExt->drhoExtFilename.resize(1);
		pl.get(e.vptInfo.drhoExt->drhoExtFilename[0], string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{	for(const string& filename: e.vptInfo.drhoExt->drhoExtFilename)
			logPrintf("%s ", filename.c_str());
	}
}
commandPerturbRhoExternal;


struct CommandPerturbIon : public Command
{
	CommandPerturbIon() : Command("perturb-ion", "jdftx/Perturb")
	{
		format = "<species> <atom> <dx0> <dx1> <dx2>";
		comments =
			"Perturb atom position. Coordinate system determined by coords-type\n";
		
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	int sp, at;
		vector3<> dx;
		pl.get(sp, 0, "species", true);
		pl.get(at, 0, "atom", true);
		pl.get(dx[0], 0.0, "dx0", true);
		pl.get(dx[1], 0.0, "dx1", true);
		pl.get(dx[2], 0.0, "dx2", true);
		
		if(e.iInfo.coordsType == CoordsLattice) dx = e.gInfo.R*dx;
		
		e.vptInfo.datom = std::make_shared<AtomPerturbation>(sp, at, dx, e);
		e.vptInfo.datom->mode.dirLattice = inv(e.gInfo.R)*e.vptInfo.datom->mode.dirCartesian; //Constructor does not initialize dirLattice properly because gInfo not set up yet
	}

	void printStatus(Everything& e, int iRep)
	{	
		auto pert = e.vptInfo.datom;
		logPrintf("%i %i %g %g %g", pert->mode.sp, pert->mode.at, pert->mode.dirCartesian[0], pert->mode.dirCartesian[1], pert->mode.dirCartesian[2]);
	}
}
commandPerturbIon;

struct CommandPerturbElectricField : public Command
{
	CommandPerturbElectricField() : Command("perturb-electric-field", "jdftx/Perturb")
	{
		format = "<Ex> <Ey> <Ez>";
		comments =
			"Perturb external electric field. This field has the same properties as the one created by the command electric-field.";
	}

	void process(ParamList& pl, Everything& e)
	{
		e.vptInfo.dElectricField = std::make_shared<ElectricFieldPerturbation>(e);
		pl.get(e.vptInfo.dElectricField->Efield[0], 0., "Ex", true);
		pl.get(e.vptInfo.dElectricField->Efield[1], 0., "Ey", true);
		pl.get(e.vptInfo.dElectricField->Efield[2], 0., "Ez", true);
	}

	void printStatus(Everything& e, int iRep)
	{	for(int k=0; k<3; k++) logPrintf("%lg ", e.vptInfo.dElectricField->Efield[k]);
	}
}
commandPerturbElectricField;


struct CommandPerturbWavevector : public Command
{
	CommandPerturbWavevector() : Command("perturb-wavevector", "jdftx/Perturb")
	{
		format = "<q0> <q1> <q2>";
		comments =
			"Set the bloch phase of the perturbation potential to exp(iGq*r) where q = (<q0>, <q1>, <q2>) is a reciprocal lattice vector\n."
			"If perturbation is not specified, will perform band structure calculation with k+q and k-q vectors.";
	}

    void process(ParamList& pl, Everything& e)
    {
		pl.get(e.vptInfo.qvec[0], 0.0, "q0");
		pl.get(e.vptInfo.qvec[1], 0.0, "q1");
		pl.get(e.vptInfo.qvec[2], 0.0, "q2");

		e.vptInfo.commensurate = !e.vptInfo.qvec.isNonzero();
	}

	void printStatus(Everything& e, int iRep)
	{
		logPrintf("%f %f %f", e.vptInfo.qvec[0], e.vptInfo.qvec[1], e.vptInfo.qvec[2]);
	}
}
commandPerturbWavevector;

struct CommandPerturbIncommensurateWavefunctions : public Command
{
	CommandPerturbIncommensurateWavefunctions() : Command("perturb-incommensurate-wavefunctions", "jdftx/Perturb")
	{
		format = "<filename> [<EcutOld>]";
		comments =
			"Incommensurate perturbations (nonzero phase as specified by set-pert-phase) require specification\n"
			"of ground state wavefunctions at wavevectors k+q and k-q for all k-points k and perturbation phase q.\n";
	}

    void process(ParamList& pl, Everything& e)
    {
		pl.get(e.vptInfo.wfnsFilename, string(), "filename", true);
		auto conversion = std::make_shared<ElecInfo::ColumnBundleReadConversion>();
		conversion->realSpace = false;
		pl.get(conversion->EcutOld, 0.0, "EcutOld");
		if(conversion->EcutOld > 0.0) e.vptInfo.readConversion = conversion;
	}

	void printStatus(Everything& e, int iRep)
	{
		if (!e.vptInfo.readConversion)
			logPrintf("%s", e.vptInfo.wfnsFilename.c_str());
		else
			logPrintf("%s %g", e.vptInfo.wfnsFilename.c_str(), e.vptInfo.readConversion->EcutOld);
	}
}
commandPerturbIncommensurateWavefunctions;

struct CommandPerturbTest : public Command
{
	CommandPerturbTest() : Command("perturb-test", "jdftx/Perturb")
	{	comments = "Perform finite difference tests of perturbation operators.\n";
	}
	void process(ParamList& pl, Everything& e) { e.vptInfo.testing = true; }
	void printStatus(Everything& e, int iRep) { return; }
}
commandPerturbTest;

