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

Energies::Energies() : TS(0), muN(0), Eband(0)
{
}


void Energies::print(FILE* fp) const
{	if(Eband)
		fprintf(fp, "Eband = %25.16lf\n\n", Eband);
	else
	{	E.print(fp, true, "%9s = %25.16lf\n");
		fprintf(fp, "-------------------------------------\n");
		fprintf(fp, "     Etot = %25.16lf\n", double(E));
		if(TS)
		{	fprintf(fp, "       TS = %25.16lf\n", TS);
			fprintf(fp, "-------------------------------------\n");
			fprintf(fp, "        F = %25.16lf\n", F());
		}
		if(muN)
		{	fprintf(fp, "      muN = %25.16lf\n", muN);
			fprintf(fp, "-------------------------------------\n");
			fprintf(fp, "        G = %25.16lf\n", G());
		}
	}
	fflush(fp);
}


double relevantFreeEnergy(const Everything& e)
{	if(e.cntrl.fixed_H) return e.ener.Eband;
	else if(e.eInfo.fillingsUpdate==ElecInfo::ConstantFillings) return double(e.ener.E);
	else if(std::isnan(e.eInfo.mu)) return e.ener.F();
	else return e.ener.G();
}
const char* relevantFreeEnergyName(const Everything& e)
{	if(e.cntrl.fixed_H) return "Eband";
	else if(e.eInfo.fillingsUpdate==ElecInfo::ConstantFillings) return "Etot";
	else if(std::isnan(e.eInfo.mu)) return "F";
	else return "G";
}



void print_Hsub_eigs(const Everything &e)
{
	const ElecInfo &eInfo = e.eInfo;
	const ElecVars &eVars = e.eVars;

	logPrintf("Band energies:\n"); //TODO: Head process should print all quantum numbers
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	logPrintf("\nstate = %d   q_k = [ %lg %lg %lg ]   w = %lg",
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
