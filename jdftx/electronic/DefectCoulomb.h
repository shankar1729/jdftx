/* Defect-only mixed-boundary electrostatics. Licensed under GPL v3 or later. */
#ifndef JDFTX_ELECTRONIC_DEFECTCOULOMB_H
#define JDFTX_ELECTRONIC_DEFECTCOULOMB_H

#include <core/Coulomb.h>

class Everything;

//! Fixed clean-system reference and variational (Kisolated-Kslab) correction.
class DefectCoulomb
{
public:
	string referenceFilename;
	double chargeTolerance = 1e-5; //!< absolute charge tolerance in electrons
	ScalarFieldTilde deltaRho; //!< unsmoothed total defect-minus-clean charge, projected to the real grid
	ScalarFieldTilde effectiveSource; //!< electron difference plus range-separated ionic difference
	ScalarFieldTilde potential, ionicPotential; //!< electron and bare-ion derivatives of the correction energy
	//! Optional isolated-kernel center in lattice coordinates; defaults to the slab embedding center.
	bool hasCenter = false;
	vector3<> center;
	double charge = 0., energy = 0.;

	void setup(const Everything& e, bool fromFixedDensity=false);
	double update(const ScalarFieldTilde& nTilde);
	void report() const;
	static void saveReference(const Everything& e, const string& filename);

private:
	const Everything* e = nullptr;
	ScalarFieldTilde nReference, ionReference;
	// Coulomb retains a reference to its parameters: these must outlive the operator.
	CoulombParams isolatedParams;
	std::shared_ptr<Coulomb> isolated;
};

#endif
