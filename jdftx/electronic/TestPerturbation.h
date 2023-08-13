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
#include <electronic/SpringConstant.h>
class Everything;
class ElecInfo;
class ColumnBundle;
class ElecVars;
class PerturbationSolver;

class TestPerturbation {
public:
	TestPerturbation(Everything &e, PerturbationSolver& ps);
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	PerturbationInfo& pInfo;
	PerturbationSolver& ps;
	SpringConstant spring;

	std::vector<ColumnBundle> Cmin;
	std::vector<ColumnBundle> Y1;
	std::vector<ColumnBundle> Y2;
	std::vector<ColumnBundle> dY;
	std::vector<ColumnBundle> C1;
	std::vector<ColumnBundle> C2;
	std::vector<ColumnBundle> C2atom;
	std::vector<ColumnBundle> dC;
	std::shared_ptr<AtomPerturbation> mode;
	std::shared_ptr<AtomPerturbation> modeA;
	std::shared_ptr<AtomPerturbation> modeB;
	std::shared_ptr<VextPerturbation> VextMode;
	vector3<> pos1;
	vector3<> pos2;

	double h = 1e-7;

	void setup(std::vector<ColumnBundle> &Y);
	void setState(std::vector<ColumnBundle> &C);
	void setAtpos1();
	void setAtpos2();
	
	void setdsqpos0();
	void setdsqposnn();
	void setdsqposnp();
	void setdsqpospn();
	void setdsqpospp();

	void testVPT();

	bool compareHamiltonians();
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
	bool FDTest_dsqV();
	bool FDTest_dsqEnl();
	bool FDTest_dsqEloc();
	bool FDTest_dsqEH();
	bool FDTest_dsqExc();
	bool FDTest_dsqExccore();
};

#endif /* ELECTRONIC_TESTPERTURBATION_H_ */
