/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <electronic/Energies.h>
#include <electronic/Everything.h>
#include <electronic/matrix.h>

Energies::Energies()
{	memset(this, 0, sizeof(*this)); 
}


void Energies::updateTotals()
{	Etot = KE + Enl + Eloc + EH + Eewald + A_diel + Eexternal + EXX + Exc + Exc_core + Epulay + EvdW;
	F = Etot - TS;
	G = F - muN;
}

void Energies::print(FILE* fp) const
{	if(Eband)
		fprintf(fp, "Eband = %25.16lf\n\n", Eband);
	else
	{	if(KE)        fprintf(fp, "KE        = %25.16lf\n", KE);
		if(Enl)       fprintf(fp, "Enl       = %25.16lf\n", Enl);
		if(Eloc)      fprintf(fp, "Eloc      = %25.16lf\n", Eloc);
		if(EH)        fprintf(fp, "EH        = %25.16lf\n", EH);
		if(Eewald)    fprintf(fp, "Eewald    = %25.16lf\n", Eewald);
		if(A_diel)    fprintf(fp, "A_diel    = %25.16lf\n", A_diel);
		if(Eexternal) fprintf(fp, "Eexternal = %25.16lf\n", Eexternal);
		if(EXX)       fprintf(fp, "EXX       = %25.16lf\n", EXX);
		if(Exc)       fprintf(fp, "Exc       = %25.16lf\n", Exc);
		if(Exc_core)  fprintf(fp, "Exc_core  = %25.16lf\n", Exc_core);
		if(Epulay)    fprintf(fp, "Epulay    = %25.16lf\n", Epulay);
		if(EvdW)      fprintf(fp, "EvdW      = %25.16lf\n", EvdW);
		fprintf(fp, "-------------------------------------\n");
		fprintf(fp, "Etot      = %25.16lf\n", Etot);
		if(TS)
		{	fprintf(fp, "TS        = %25.16lf\n", TS);
			fprintf(fp, "-------------------------------------\n");
			fprintf(fp, "F         = %25.16lf\n", F);
		}
		if(muN)
		{	fprintf(fp, "muN       = %25.16lf\n", muN);
			fprintf(fp, "-------------------------------------\n");
			fprintf(fp, "G         = %25.16lf\n", G);
		}
	}
	fflush(fp);
}


double relevantFreeEnergy(const Everything& e)
{	if(e.cntrl.fixed_n) return e.ener.Eband;
	else if(e.eInfo.fillingsUpdate==ElecInfo::ConstantFillings) return e.ener.Etot;
	else if(isnan(e.eInfo.mu)) return e.ener.F;
	else return e.ener.G;
}
const char* relevantFreeEnergyName(const Everything& e)
{	if(e.cntrl.fixed_n) return "Eband";
	else if(e.eInfo.fillingsUpdate==ElecInfo::ConstantFillings) return "Etot";
	else if(isnan(e.eInfo.mu)) return "F";
	else return "G";
}



void print_Hsub_eigs(const Everything &e)
{
	const ElecInfo &eInfo = e.eInfo;
	const ElecVars &eVars = e.eVars;

	logPrintf("Band energies:\n");
	for(int q=0; q < eInfo.nStates; q++)
	{
		logPrintf("\nstate = %d   q_k = [ %lg %lg %lg ]   w = %lg",
			q, eInfo.qnums[q].k[0],eInfo.qnums[q].k[1], eInfo.qnums[q].k[2], eInfo.qnums[q].weight);
		logPrintf("   spin = %d\n",eInfo.qnums[q].spin);
		logPrintf("%4s  %13s  %13s  %13s\n",
			"band","filling   ","diag(Hsub) ","epsilon   ");
		logPrintf("-------------------------------------------------\n");
		//NOTE: fillings are always 0 to 1 internally, but read/write 0 to 2 for SpinNone
		double spinFactor =  (eInfo.spinType==SpinNone ? 2 : 1);
		const diagMatrix diagHsubq = diag(eVars.Hsub[q]);
		for(int i=0; i < eInfo.nBands; i++)
			logPrintf("%4d  %13.6le  %13.6le  %13.6le\n",
				i, spinFactor*eVars.F[q][i], diagHsubq[i], eVars.Hsub_eigs[q][i]);
	}
	logFlush();
}

