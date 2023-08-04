/*
 * PerturbationSolver.h
 *
 *  Created on: Jul 22, 2022
 *      Author: brandon
 */

#ifndef ELECTRONIC_PERTURBATIONSOLVER_H_
#define ELECTRONIC_PERTURBATIONSOLVER_H_


#include <core/Minimize.h>
#include <core/ScalarFieldArray.h>

class Everything;
class ElecInfo;
class ColumnBundle;
class ElecVars;
class PerturbationInfo;
class Basis;
class QuantumNumber;
class IonInfo;

struct PerturbationGradient {
	std::vector<ColumnBundle> dY; //!< wavefunctions
	const ElecInfo* eInfo;
	const PerturbationInfo* pInfo;

	void init(const Everything& e); //!< initialize C and Haux with the correct sizes for everything
	PerturbationGradient& operator*=(double alpha);
	//PerturbationGradient& operator+=(const IonicGradient&);
	//PerturbationGradient& operator-=(const IonicGradient&);

	//PerturbationGradient operator*(double alpha) const;
	//PerturbationGradient operator+(const IonicGradient&) const;
	//PerturbationGradient operator-(const IonicGradient&) const;
};

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y); //!< accumulate operation: Y += alpha*X
double dot(const PerturbationGradient& x, const PerturbationGradient& y); //!< inner product
double dot(const std::vector<ColumnBundle>& x, const std::vector<ColumnBundle>& y, PerturbationInfo* pInfo, ElecInfo *eInfo);
PerturbationGradient clone(const PerturbationGradient& x); //!< create a copy
void randomize(PerturbationGradient& x); //!< initialize with random numbers
ScalarField randomRealSpaceVector(ColumnBundle& ref, const  GridInfo* gInfo);
void printVvpt(ScalarField V);

class PerturbationSolver : public Minimizable<PerturbationGradient>
{
public:
	PerturbationSolver(Everything& e);

	void randomize(PerturbationGradient& x);
	void printCB(ColumnBundle C);
	void printM(matrix C);
	void printV(ScalarField V);
	void printV(ScalarFieldTilde V);

	double compute(PerturbationGradient* grad, PerturbationGradient* Kgrad);
	void getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y);
	void step(const PerturbationGradient& dir, double alpha);
	void constrain(PerturbationGradient&);
	double minimize(const MinimizeParams& params);
	void computeIncommensurateWfns();

	ScalarFieldArray getdn(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C = 0);
	ScalarFieldArray getdnatom();
	void getdnInc(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq);
	ScalarFieldArray getn(const std::vector<ColumnBundle>& C);

	void updateHC();
	void updateNonlocalDerivs();
	void applyH(const QuantumNumber& qnum, const diagMatrix& Fq, ColumnBundle& HCq, const ColumnBundle& Cq);
	//void dAugatom(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C);
	void dH_Vscloc(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq = 0);
	void dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq = 0);
	void dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const complexScalarFieldArray& dVscloc);
	void dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc);
	void dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, int qsign);
	ScalarField get_nTot(const ScalarFieldArray n);
	ScalarFieldArray getVscloc(ScalarFieldArray ntot);
	void getdVsclocPsi(const ScalarFieldArray dn, ScalarFieldArray& dVscloc);
	void getdVsclocPsi(const complexScalarFieldArray dn, complexScalarFieldArray& dVscloc, vector3<>* q);
	void getdVsclocTau(ScalarFieldArray& dVscloc);
	void getdVsclocTau(complexScalarFieldArray& dVscloc, int qsign);
	void orthonormalize(int q, ColumnBundle& Cq);
	void calcdGradTau();
	ScalarFieldArray addnXC(const ScalarFieldArray n);
	ScalarFieldArray getdnXC(const ScalarField dnCore);
    
    bool ultrasoftDensityAugRequired = false;
	
	ScalarFieldTilde dVlocpsA;
	ScalarFieldArray dVxcA;
	ScalarField dnCoreA;
private:
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	IonInfo& iInfo;
	PerturbationInfo& pInfo;
};

void init(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const Basis* basis, const ElecInfo* eInfo);

#endif /* ELECTRONIC_PERTURBATIONSOLVER_H_ */
