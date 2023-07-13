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

	std::vector<string> dVexternalFilename; //!< external potential filename (read in real space)
	string wfnsFilename;

	std::vector<ColumnBundle> dY, dC, dGradTau, dGradPsi;
	std::vector<matrix> dU;
	complexScalarFieldArray dVext;
	ScalarFieldArray dn;
	complexScalarFieldArray dnpq;
	complexScalarFieldArray dnmq;

	std::vector<ColumnBundle> Cinc, HC;

	vector3<> qvec;
	bool incommensurate = false;
	std::vector<QuantumNumber> Tk_vectors;
	std::vector<QuantumNumber> Tinvk_vectors;
	std::vector<Basis> Tk_basis;
	std::vector<Basis> Tinvk_basis;

	bool testing = false;
};

#endif /* ELECTRONIC_PERTURBATIONINFO_H_ */

























