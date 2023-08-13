/*
 * PerturbationInfo.h
 *
 *  Created on: Jul 25, 2022
 *      Author: brandon
 */

#ifndef ELECTRONIC_PERTURBATIONINFO_H_
#define ELECTRONIC_PERTURBATIONINFO_H_

#include <core/ScalarField.h>
#include <electronic/ColumnBundle.h>
#include <vector>
#include <string>
#include <core/matrix.h>
#include <electronic/ElecInfo.h>
#include <electronic/ExCorr.h>
#include <electronic/PerturbationSolver.h>

class Everything;
class ElecVars;

struct AtomicMode
{	int sp, at; //!< species and atom number (within first unit cell)
	vector3<> dirLattice; //!< dx in lattice coords
	vector3<> dirCartesian; //!< dx in Cartesian coords
};

class AtomPerturbation {
public:
	AtomPerturbation(int sp, int at, int iDir, const Everything& e);
	AtomPerturbation(int sp, int at, vector3<> dirCartesian, const Everything& e);
	bool sameAtom(const std::shared_ptr<AtomPerturbation> pert);
	void init(const Everything &e, const ElecVars& eVars, const ElecInfo& eInfo);
    bool isUltrasoft(const IonInfo& iInfo); //!< Does species use ultrasoft pp
	
	AtomicMode mode; //Used for atomic perturbations
	std::vector<ColumnBundle> Vatom_cached, dVatom_cached;
	std::vector<matrix> VdagCatom_cached, dVdagCatom_cached;
	matrix E_nAug_datom;
	ScalarFieldTilde Vlocps;
	ScalarFieldArray dnatom;
	std::vector<ColumnBundle> dCatom;
};

class VextPerturbation {
public:
	std::vector<string> dVexternalFilename; //!< External potential filename (read in real space)
	ScalarFieldArray dVext;
	complexScalarFieldArray dVextpq, dVextmq;
};

class PerturbationInfo {
public:
	PerturbationInfo();
	virtual ~PerturbationInfo();
	void setup(const Everything &e, const ElecVars &eVars);
	void read(const ElecInfo &eInfo, const Everything &e, std::vector<ColumnBundle>& C, const char *fname, const ElecInfo::ColumnBundleReadConversion* conversion=0) const;
	void setupkpoints(const Everything &e, const ElecInfo &eInfo);
	void initInc(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const ElecInfo* eInfo); //!< Setup incommensurate ColumnBundles
	void checkSupportedFeatures(const Everything &e, const ElecInfo &eInfo);
	bool densityAugRequired(const Everything &e);

	bool commensurate = true; //!< Is the perturbation lattice periodic
	bool testing = false; //!< Whether or not a FD test is being conducted
	
	std::shared_ptr<VextPerturbation> dVext; //!< Is there a perturbation of Vext
    std::shared_ptr<AtomPerturbation> datom;
	
	vector3<> qvec; //!< Bloch wavevector of perturbation. Equal to zero if perturbation is commensurate
	
	std::vector<QuantumNumber> Tk_vectors; //!< List of k+q vectors
	std::vector<QuantumNumber> Tinvk_vectors; //!< List of k-q vectors
	std::vector<Basis> Tk_basis; //!< List of k+q basis elements
	std::vector<Basis> Tinvk_basis; //!< List of k-q basis elements
	
	string wfnsFilename; //!< Name of ground state wavefunctions

	/* Gradients and wavefunction shifts */
	PerturbationGradient dGradTau, dGradPsi;
	std::vector<ColumnBundle> dC, Cinc;
	std::vector<matrix> dU, dUmhalfatom, dHsub, dHsubatom;
    
    /* Intermediate scalar field derivs */
	ScalarFieldArray dn;
	ScalarField dnCoreA;
	complexScalarFieldArray dnpq, dnmq;
	ScalarFieldArray dVscloc, dVsclocTau;
	complexScalarFieldArray dVsclocpq, dVsclocmq;
    
    /* Cached quantities */
	std::vector<ColumnBundle> grad, HC, OC;
    ScalarField sigma_cached, e_nn_cached, e_sigma_cached, e_nsigma_cached, e_sigmasigma_cached; //LDA and GGA second derivs
    VectorField IDJn_cached;
	
	std::vector<matrix> E_nAug_cached, E_nAug_dVsclocpsi, E_nAug_dVscloctau;
};

#endif /* ELECTRONIC_PERTURBATIONINFO_H_ */

























