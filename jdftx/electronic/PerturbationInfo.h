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

class Everything;
class ElecInfo;
class ElecVars;
//class ColumnBundle;
//#include <electronic/Everything.h>

class PerturbationInfo {
public:
	PerturbationInfo();
	virtual ~PerturbationInfo();
	void setup(Everything *e, ElecVars *eVars);

	std::vector<string> dVexternalFilename; //!< external potential filename (read in real space)
	std::vector<ColumnBundle> dY, dC, dGradTau, dGradPsi;
	std::vector<matrix> dU;
	std::vector<ScalarField> dVext;
	//TODO q-vector
};

#endif /* ELECTRONIC_PERTURBATIONINFO_H_ */

























