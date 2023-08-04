/*
 * TestPerturbation.h
 *
 *  Created on: Oct 30, 2022
 *      Author: brandon
 */

#ifndef ELECTRONIC_TESTPERTURBATION_H_
#define ELECTRONIC_TESTPERTURBATION_H_

#include <vector>
#include <electronic/PerturbationInfo.h>
class Everything;
class ElecInfo;
class ColumnBundle;
class ElecVars;
class PerturbationSolver;

class TestPerturbation {
public:
	TestPerturbation(Everything &e);
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	PerturbationInfo& pInfo;
	PerturbationSolver* ps = 0;

	std::vector<ColumnBundle> Cmin;
	std::vector<ColumnBundle> Y1;
	std::vector<ColumnBundle> Y2;
	std::vector<ColumnBundle> dY;
	std::vector<ColumnBundle> C1;
	std::vector<ColumnBundle> C2;
	std::vector<ColumnBundle> C2atom;
	std::vector<ColumnBundle> dC;
	PerturbationInfo::Mode m;
	vector3<> pos1;
	vector3<> pos2;

	double h = 1e-7;

	void setup(std::vector<ColumnBundle> &Y);
	void setState(std::vector<ColumnBundle> &C);
	void setAtpos1();
	void setAtpos2();

	void testVPT();

	bool compareHamiltonians();
	bool compareVxc();
	bool compareVscloc();
	bool compareGetn();
	bool testGradientIsZero();
	bool FDTest_dVxc();
	bool FDTest_dn();
	bool FDTest_dVscloc();
	bool FDTest_dnatom();
	bool FDTest_dVsclocatom();
	bool FDTest_dVlocpsatom();
	bool FDTest_dgradpsi();
	bool FDTest_Hamiltonian();
	bool FDTest_Hamiltoniandatom();
	bool FDTest_Overlapatom();
	bool FDTest_dV();
	bool FDTest_dC();
	bool FDTest_dCatom();
	bool FDTest_dgradtau();
	bool FDTest_dgradtauatom();
};

#endif /* ELECTRONIC_TESTPERTURBATION_H_ */
