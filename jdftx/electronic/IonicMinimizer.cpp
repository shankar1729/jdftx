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

#include <electronic/IonicMinimizer.h>
#include <electronic/IonInfo.h>
#include <electronic/Symmetries.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/Dump.h>
#include <electronic/operators.h>
#include <core/Random.h>
#include <core/BlasExtra.h>

void IonicGradient::init(const IonInfo& iInfo)
{	clear();
	resize(iInfo.species.size());
	for(unsigned sp=0; sp<size(); sp++)
		at(sp).resize(iInfo.species[sp]->atpos.size());
}

void IonicGradient::print(const Everything& e, FILE* fp, const char* prefix) const
{	fprintf(fp, "# Forces in %s coordinates:\n", forcesOutputCoordsMap.getString(e.iInfo.forcesOutputCoords));
	for(unsigned sp=0; sp<size(); sp++)
	{	const SpeciesInfo& sinfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<at(sp).size(); atom++)
		{	vector3<> ff;
			switch(e.iInfo.forcesOutputCoords)
			{	case ForcesCoordsLattice: ff = at(sp)[atom]; break;
				case ForcesCoordsCartesian: ff = e.gInfo.invRT * at(sp)[atom]; break;
				case ForcesCoordsContravariant: ff = e.gInfo.invRTR * at(sp)[atom]; break;
				case ForcesCoordsPositions: assert(false); //should not get here
			}
			fprintf(fp, "%s %s %19.15lf %19.15lf %19.15lf %lg", prefix,
				sinfo.name.c_str(), ff[0], ff[1], ff[2], sinfo.constraints[atom].moveScale);
			if(sinfo.constraints[atom].type != SpeciesInfo::Constraint::None)
				sinfo.constraints[atom].print(fp, e);
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n");
}


IonicGradient& IonicGradient::operator*=(double s)
{	for(unsigned sp=0; sp<size(); sp++)
		for(unsigned atom=0; atom<at(sp).size(); atom++)
			at(sp)[atom] *= s;
	return *this;
}

IonicGradient IonicGradient::operator*(double s) const
{	IonicGradient result(*this);
	result *= s;
	return result;
}

IonicGradient& IonicGradient::operator+=(const IonicGradient& other)
{	axpy(1., other, *this);
	return *this;
}

IonicGradient IonicGradient::operator+(const IonicGradient& other) const
{	IonicGradient result(*this);
	axpy(1., other, result);
	return result;
}

IonicGradient IonicGradient::operator-(const IonicGradient& other) const
{	IonicGradient result(*this);
	axpy(-1., other, result);
	return result;
}

void axpy(double alpha, const IonicGradient& x, IonicGradient& y)
{	assert(x.size() == y.size());
	for(unsigned sp=0; sp<x.size(); sp++)
	{	assert(x[sp].size() == y[sp].size());
		for(unsigned atom=0; atom<x[sp].size(); atom++)
			y[sp][atom] += alpha * x[sp][atom];
	}
}

double dot(const IonicGradient& x, const IonicGradient& y)
{	double result = 0.0;
	assert(x.size() == y.size());
	for(unsigned sp=0; sp<x.size(); sp++)
	{	assert(x[sp].size() == y[sp].size());
		for(unsigned atom=0; atom<x[sp].size(); atom++)
			result += dot(x[sp][atom], y[sp][atom]);
	}
	mpiUtil->bcast(result); //ensure consistency to numerical precision
	return result;
}

IonicGradient clone(const IonicGradient& x)
{	return x; //implicit copy constructor handles everything correctly
}

void randomize(IonicGradient& x)
{	for(unsigned sp=0; sp<x.size(); sp++)
		for(unsigned atom=0; atom<x[sp].size(); atom++)
			for(int k=0; k<3; k++)
				x[sp][atom][k] = Random::normal();
}


IonicGradient operator*(const matrix3<>& mat, const IonicGradient& x)
{	IonicGradient ret(x);
	for(unsigned sp=0; sp<x.size(); sp++)
		for(unsigned atom=0; atom<x[sp].size(); atom++)
			ret[sp][atom] = mat * x[sp][atom];
	return ret;
}


IonicMinimizer::IonicMinimizer(Everything& e) : e(e)
{	
}

void IonicMinimizer::step(const IonicGradient& dir, double alpha)
{	static StopWatch watch("WavefunctionDrag"); watch.start();
	ElecVars& eVars = e.eVars;
	IonInfo& iInfo = e.iInfo;
	
	IonicGradient dpos = alpha * e.gInfo.invR * dir; //dir is in cartesian, atpos in lattice
	
	if(e.cntrl.dragWavefunctions)
	{	//Check if atomic orbitals available and compile list of displacements for each orbital:
		std::vector< vector3<> > drColumns;
		for(unsigned s=0; s<iInfo.species.size(); s++)
		{	const SpeciesInfo& sp = *(iInfo.species[s]);
			int nOrb = sp.nAtomicOrbitals() / sp.atpos.size(); //number of orbitals per atom
			for(int iOrb=0; iOrb<nOrb; iOrb++)
				drColumns.insert(drColumns.end(), dpos[s].begin(), dpos[s].end());
		}
		
		if(drColumns.size()) 
		{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			{
				//Get atomic orbitals at old positions:
				eVars.Y[q].free();
				ColumnBundle psi = iInfo.getAtomicOrbitals(q);
				
				//Fit the wavefunctions to atomic orbitals (minimize C0^OC0 where C0 is the remainder)
				matrix alpha;  //LCAO coefficients for best fit
				{	ColumnBundle Opsi = O(psi); //non-trivial cost for uspp
					alpha = inv(psi^Opsi) * (Opsi^eVars.C[q]);
				}
				eVars.C[q] -= psi * alpha; //now contains residual
				
				//Translate the atomic orbitals and reconsitute wavefunctions:
				translateColumns(psi, drColumns.data());
				eVars.C[q] += psi * alpha;
				
				//Orthonormalize and replace old wavefunctions:
				eVars.Y[q] = eVars.C[q] * invsqrt(eVars.C[q]^O(eVars.C[q]));
			}
		}
	}

	//Move the atoms:
	for(unsigned sp=0; sp < iInfo.species.size(); sp++)
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.atpos.size(); atom++)
			spInfo.atpos[atom] += dpos[sp][atom]; 
		mpiUtil->bcast((double*)spInfo.atpos.data(), 3*spInfo.atpos.size());
		#ifdef GPU_ENABLED
		spInfo.sync_atposGpu();
		#endif
	}
	watch.stop();
}

double IonicMinimizer::compute(IonicGradient* grad)
{
	//Initialize ion-dependent quantities at this position:
	e.iInfo.update(e.ener);

	//Minimize the electronic system:
	elecFluidMinimize(e);
	
	//Calculate forces if needed:
	if(grad)
	{	e.iInfo.ionicEnergyAndGrad(e.iInfo.forces); //compute forces in lattice coordinates
		*grad = -e.gInfo.invRT * e.iInfo.forces; //gradient in cartesian coordinates (and negative of force)
	}
	
	double ener = relevantFreeEnergy(e);
	mpiUtil->bcast(ener);
	return ener;
}

IonicGradient IonicMinimizer::precondition(const IonicGradient& grad)
{	IonicGradient Kgrad(grad);
	for(unsigned sp=0; sp<grad.size(); sp++)
	{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<grad[sp].size(); atom++)
			Kgrad[sp][atom] = spInfo.constraints[atom](grad[sp][atom]);  //Apply move constraints to the force gradient.
	}
	return Kgrad;
}

bool IonicMinimizer::report(int iter)
{	e.iInfo.printPositions(globalLog);
	e.iInfo.forces.print(e, globalLog);
	e.dump(DumpFreq_Ionic, iter);
	return false;
}

void IonicMinimizer::constrain(IonicGradient& x)
{	e.symm.symmetrize(x); //minimization is constrained by symmetries on forces
}
