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

	double compute(PerturbationGradient* grad, PerturbationGradient* Kgrad);
	void getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y);
	void step(const PerturbationGradient& dir, double alpha);
	bool report(int iter);
	void constrain(PerturbationGradient&);
	double minimize(const MinimizeParams& params);
	void computeIncommensurateWfns();

	ScalarFieldArray getdn(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C = 0);
	void getdnInc(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq);
	ScalarFieldArray getn(std::vector<ColumnBundle>& C);

	void applyH(QuantumNumber qnum, const diagMatrix& F, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray ntot);
	void applyH2(int q, const diagMatrix& Fq, ColumnBundle& HCq);
	void dH_Vscloc(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dVscloc);
	void dHpsi(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dn);
	void dHpsi(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, complexScalarFieldArray dn, vector3<> q);
	void dHtau(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C);
	void dHtau(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, int qsign);
	ScalarField get_nTot(ScalarFieldArray n);
	ScalarFieldArray get_nXC(ScalarFieldArray n);
	ScalarFieldArray getVscloc(ScalarFieldArray ntot);
	void getdVsclocPsi(ScalarFieldArray dn, ScalarFieldArray& dVscloc);
	void getdVsclocPsi(complexScalarFieldArray dn, complexScalarFieldArray& dVscloc, vector3<> k);
	void getdVsclocTau(complexScalarFieldArray& dVscloc, bool conjugate=false);
	void orthonormalize(int q, ColumnBundle& Cq);
	void calcdGradTau();
	ScalarFieldArray addnXC(ScalarFieldArray n);
private:
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	PerturbationInfo& pInfo;
};

void init(std::vector<ColumnBundle>& Y, int nbundles, int ncols, const Basis* basis, const ElecInfo* eInfo);

#endif /* ELECTRONIC_PERTURBATIONSOLVER_H_ */
