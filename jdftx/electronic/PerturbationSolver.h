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

struct PerturbationGradient {
	std::vector<ColumnBundle> dY; //!< wavefunctions
	const ElecInfo* eInfo;

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


class PerturbationSolver : public Minimizable<PerturbationGradient>
{
public:
	PerturbationSolver(Everything& e);

	double compute(PerturbationGradient* grad, PerturbationGradient* Kgrad);
	void getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y);
	void step(const PerturbationGradient& dir, double alpha);
	bool report(int iter);
	void constrain(PerturbationGradient&);
	double minimize(const MinimizeParams& params);

	ScalarFieldArray getdn(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C = 0);
	ScalarFieldArray getn(std::vector<ColumnBundle>& C);

	void applyH(int q, const diagMatrix& F, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray ntot);
	void applyH2(int q, const diagMatrix& Fq, ColumnBundle& HCq);
	void dHpsi(ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dn);
	void dHtau(ColumnBundle& HC, ColumnBundle& C);
	ScalarField get_nTot(ScalarFieldArray n);
	ScalarFieldArray get_nXC(ScalarFieldArray n);
	ScalarFieldArray getVscloc(ScalarFieldArray ntot);
	void getdVsclocPsi(ScalarFieldArray dn, ScalarFieldArray& dVscloc);
	void getdVsclocTau(ScalarFieldArray& dVscloc);
	void orthonormalize(int q, ColumnBundle& Cq);
	void calcdGradTau();
	void testVPT();
	ScalarFieldArray addnXC(ScalarFieldArray n);
private:
	Everything& e;
	ElecVars& eVars;
	const ElecInfo& eInfo;
	PerturbationInfo& pInfo;
};


#endif /* ELECTRONIC_PERTURBATIONSOLVER_H_ */
