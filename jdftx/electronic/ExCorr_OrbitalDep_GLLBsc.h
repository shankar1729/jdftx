/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_EXCORR_ORBITALDEP_GLLBSC_H
#define JDFTX_ELECTRONIC_EXCORR_ORBITALDEP_GLLBSC_H

#include <electronic/ExCorr.h>

struct ExCorr_OrbitalDep_GLLBsc : public ExCorr::OrbitalDep
{	ExCorr_OrbitalDep_GLLBsc(const Everything&);
	bool ignore_nCore() const { return true; }
	ScalarFieldArray getPotential() const;
	void dump() const;
private:
	double T; //smearing width
	std::vector<double> getExtremalEnergy(bool HOMO) const; //!<  get HOMO or LUMO energy (depending on HOMO=true/false), optionally accounting for smearing (depending on T)
	ScalarFieldArray getPotential(std::vector<double> eHOMO, std::vector<double>* eLUMO=0) const; //!< get the orbital dep potential (or discontinuity contribution if eLUMO is non-null)
};

#endif //JDFTX_ELECTRONIC_EXCORR_ORBITALDEP_GLLBSC_H
