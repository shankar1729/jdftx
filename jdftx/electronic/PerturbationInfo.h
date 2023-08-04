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

class Everything;
//class ElecInfo;
class ElecVars;
//class ColumnBundleReadConversion;

//class ColumnBundle;
//#include <electronic/Everything.h>

class PerturbationInfo {
public:
	PerturbationInfo();
	virtual ~PerturbationInfo();
	void setup(const Everything &e, const ElecVars &eVars);
	void read(const ElecInfo &eInfo, const Everything &e, std::vector<ColumnBundle>& C, const char *fname, const ElecInfo::ColumnBundleReadConversion* conversion=0) const;
	void setupkpoints(const Everything &e, const ElecInfo &eInfo);
	void initInc(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const ElecInfo* eInfo);
    void updateExcorrCache(const ExCorr& exc, const GridInfo& gInfo, const ScalarField& n);

	std::vector<string> dVexternalFilename; //!< external potential filename (read in real space)
	string wfnsFilename;

	std::vector<ColumnBundle> dGradTau, dGradPsi;
	std::vector<ColumnBundle> dY, dC, HdC, dHC, Cinc;
	std::vector<ColumnBundle> dCatom, HdCatom, dHCatom;
	std::vector<matrix> dU, dUsqrtinvatom;
    
	ScalarFieldArray dVext;
	complexScalarFieldArray dVextpq, dVextmq;
    bool VextPerturbationExists = false;
    bool atposPerturbationExists = false;
    
	ScalarFieldArray dn, dnatom;
	complexScalarFieldArray dnpq, dnmq;
	ScalarFieldArray dVscloc, dVsclocatom;
	complexScalarFieldArray dVsclocpq, dVsclocmq;
    
    /* Cached quantities */
	std::vector<ColumnBundle> HC, OC;
    ScalarField sigma_cached, e_nn_cached, e_sigma_cached, e_nsigma_cached, e_sigmasigma_cached; //Excorr second derivs
    VectorField IDJn_cached;
	
	std::vector<ColumnBundle> Vatom_cached, dVatom_cached;
	std::vector<matrix> VdagCatom_cached, dVdagCatom_cached;
	
	std::vector<matrix> E_nAug_cached, E_nAug_dVsclocpsi, E_nAug_dVscloctau;
	matrix E_nAug_datom;

	vector3<> qvec;
	std::vector<QuantumNumber> Tk_vectors;
	std::vector<QuantumNumber> Tinvk_vectors;
	std::vector<Basis> Tk_basis;
	std::vector<Basis> Tinvk_basis;

	bool incommensurate = false;
	bool testing = false;
	
	struct Mode
	{	int sp, at; //!< species and atom number (within first unit cell)
		vector3<> dirLattice; //!< 
		vector3<> dirCartesian; //!< 
	};
	
	Mode atomdisplacement;
};

#endif /* ELECTRONIC_PERTURBATIONINFO_H_ */

























