/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman

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

#ifndef JDFTX_WANNIER_DEFECT_SUPERCELL_H
#define JDFTX_WANNIER_DEFECT_SUPERCELL_H

#include <electronic/ColumnBundle.h>
#include <memory>

class DefectSupercell
{
public:
	string name; //defect name
	string inFile; //input file of supercell calculation
	vector3<int> supIn; //supercell of defect calculation relative to unit cell (determined automatically and need not be commensurate)
	vector3<int> supOut; //supercell for Wannierized matrix element output (must be commensurate with k-point supercell)
	vector3<> xCenter; //lattice coordinates of supercell where defect is centered
	double q; //net electron count of defect
	
	void initialize(const class Wannier*);
	
	struct CachedProjections
	{	std::vector<matrix> VdagC[2]; //projections for reference and defect supercell positions
		std::vector<matrix> psiDagC[2]; //corresponding atomic-orbital projections for DFT+U, if needed
	};
	void project(const ColumnBundle& C, CachedProjections& proj); //compute projections for given wavefunctions
	void bcast(CachedProjections& proj, int src) const; //make projections available on all processes
	
	//Return defect matrix elements between C1 and C2 (at different k)
	matrix compute(const ColumnBundle& C1, const ColumnBundle& C2,
		const CachedProjections& proj1, const CachedProjections& proj); //corresponding  projections

private:
	const class Wannier* wannier;
	const class Everything* e; //unit cell calculation
	std::shared_ptr<struct WannierEverything> eSup; //defect calculation
	std::shared_ptr<class WignerSeitz> wsSup; //Wigner-Seitz cell of supercell
	std::vector<std::vector<vector3<>>> atposRef; //unit cell atomic positions by eSup species in eSup Wigner-Seitz cell
	std::vector<matrix> U_rhoAtomRef; //e.eVars.U_rhoAtom repeated appropriately for perfect supercell
	std::vector<matrix> Urho[2]; //DFT+U potentials by species for reference and defect supercells
};

#endif //JDFTX_WANNIER_DEFECT_SUPERCELL_H
