/* Defect-only mixed-boundary electrostatics. Licensed under GPL v3 or later. */
#include <commands/command.h>
#include <electronic/Everything.h>
#include <electronic/DefectCoulomb.h>
#include <cmath>

struct CommandDefectCoulomb : public Command
{
	CommandDefectCoulomb() : Command("defect-coulomb", "jdftx/Coulomb interactions")
	{	format = "<referenceFile> [<chargeTolerance>=1e-5]";
		comments = "Apply (Kisolated-Kslab) to the current total charge minus a fixed clean reference\n"
			"at every electronic density/potential evaluation, with variational energy correction.\n"
			"Generate <referenceFile> using dump End DefectReference in a converged neutral clean slab.\n"
			"Requires identical cell, FFT grid, cutoffs, spin, smearing and host pseudopotentials.\n"
			"Requires electronic-scf, symmetries none, vacuum, fixed neutral electron number,\n"
			"fixed ions/lattice and Slab Coulomb interaction, with optional embedding.\n"
			"Embedded references require the same embedding center and ionic margin.\n"
			"The isolated operator inherits the slab embedding center unless defect-coulomb-center is set.\n"
			"The density difference\n"
			"must satisfy the confinement requirements of Isolated WS truncation.\n"
			"Charge tolerance is in electrons; a non-neutral difference aborts without adding a background.\n"
			"dump End DefectCoulomb writes deltaRho and Vcorr for localization diagnostics.";
		require("coulomb-interaction");
		require("electronic-scf");
	}
	void process(ParamList& pl, Everything& e)
	{	auto correction = std::make_shared<DefectCoulomb>();
		pl.get(correction->referenceFilename, string(), "referenceFile", true);
		pl.get(correction->chargeTolerance, 1e-5, "chargeTolerance");
		if(!correction->referenceFilename.length()) throw string("<referenceFile> must not be empty.");
		if(!std::isfinite(correction->chargeTolerance) || correction->chargeTolerance<=0.)
			throw string("<chargeTolerance> must be finite and positive.");
		e.eVars.defectCoulomb = correction;
	}
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s %lg", e.eVars.defectCoulomb->referenceFilename.c_str(), e.eVars.defectCoulomb->chargeTolerance);
	}
} commandDefectCoulomb;

struct CommandDefectCoulombCenter : public Command
{
	CommandDefectCoulombCenter() : Command("defect-coulomb-center", "jdftx/Coulomb interactions")
	{	format = "<c0> <c1> <c2>";
		comments = "Override the center of the isolated defect Coulomb operator only.\n"
			"The slab operator retains its coulomb-truncation-embed center.\n"
			"Coordinates follow coords-type (Cartesian values are in bohr).\n"
			"Requires embedded truncation; choose the lateral center near the localized defect.\n"
			"Default: inherit all three coordinates of coulomb-truncation-embed.";
		require("defect-coulomb"); require("coulomb-truncation-embed");
		require("coords-type"); require("latt-scale");
	}
	void process(ParamList& pl, Everything& e)
	{	auto& correction = *e.eVars.defectCoulomb;
		for(int k=0; k<3; k++)
		{	pl.get(correction.center[k], 0., "center coordinate", true);
			if(!std::isfinite(correction.center[k])) throw string("Defect center coordinates must be finite.");
		}
		if(e.iInfo.coordsType==CoordsCartesian) correction.center = inv(e.gInfo.R)*correction.center;
		correction.hasCenter = true;
	}
	void printStatus(Everything& e, int iRep)
	{	vector3<> center = e.eVars.defectCoulomb->center;
		if(e.iInfo.coordsType==CoordsCartesian) center = e.gInfo.R*center;
		logPrintf("%.15g %.15g %.15g", center[0], center[1], center[2]);
	}
} commandDefectCoulombCenter;
