/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef FLUID1D_CORE1D_DATACOLLECTION_H
#define FLUID1D_CORE1D_DATACOLLECTION_H


//! @file DataCollection.h
//! @brief classes ScalarFieldCollection, ScalarFieldTildeCollection and just enough operators to enable CG w.r.t to them

#include <core/Util.h>
#include <core1D/Operators.h>
#include <vector>

#include <cstdio>

typedef std::vector<ScalarField> ScalarFieldCollection; //!< dynamic size collection of real-space scalar fields
typedef std::vector<ScalarFieldTilde> ScalarFieldTildeCollection; //!< dynamic size collection of basis-space scalar fields
#define TCollection std::vector<T> //!< shorthand for templates below (undef'd at end of file)

//! Scale
template<typename T> TCollection& operator*=(TCollection& x, double alpha)
{ 	for(unsigned i=0; i<x.size(); i++) x[i] *= alpha;
	return x;
}

//! y += alpha x
template<typename T> void axpy(double alpha, const TCollection& x, TCollection& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) axpy(alpha, x[i], y[i]);
}

//! Inner product
template<typename T> double dot(const TCollection& x, const TCollection& y)
{	assert(x.size()==y.size());
	double ret = 0.0;
	for(unsigned i=0; i<x.size(); i++) ret += dot(x[i], y[i]);
	return ret;
}

//! Create a copy of the data (note operator= references same data since Tptr's are pointers!)
template<typename T> TCollection clone(const TCollection& x)
{	TCollection ret(x.size());
	for(unsigned i=0; i<x.size(); i++) ret[i] = x[i];
	return ret;
}

//! Initialize (non-null) data to zero
template<typename T> void initZero(TCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initZero(x[i]);
}

//! Allocate all null components and initialize them to zero
//! If N is non-zero, resize collection to N scalar fields
template<typename T> void nullToZero(TCollection& x, const GridInfo& gInfo, int N=0)
{	if(N) x.resize(N);
	for(unsigned i=0; i<x.size(); i++) nullToZero(x[i], gInfo);
}

//! Initialize to random numbers (uniform on 0 to 1)
template<typename T> void initRandomFlat(TCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandomFlat(x[i]);
}
//! Initialize to normal random numbers:
template<typename T> void randomize(TCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandom(x[i], 3.0);
}

template<typename T> void loadFromFile(TCollection& x, const char* filename)
{	FILE* fp = fopen(filename, "rb");
	if(!fp) die("Could not open %s for reading.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in loadFromFile(x,\"%s\").\n", i, filename)
		if(fread(x[i].data(), sizeof(double), x[i].nData(), fp) < unsigned(x[i].nData()))
			die("File ended too soon while reading x[%d] in loadFromFile(x,\"%s\").\n", i, filename)
	}
	fclose(fp);
}

template<typename T> void saveToFile(const TCollection& x, const char* filename)
{	FILE* fp = fopen(filename, "wb");
	if(!fp) die("Could not open %s for writing.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in saveToFile(x,\"%s\").\n", i, filename)
		fwrite(x[i].data(), sizeof(double), x[i].nData(), fp);
	}
	fclose(fp);
}

inline const std::vector<double>& getRadialCoordinate(const ScalarField& X) { return X.gInfo->r; }
inline const std::vector<double>& getRadialCoordinate(const ScalarFieldTilde& X) { return X.gInfo->G; }

template<typename T> void printToFile(const TCollection& x, const char* filename)
{	//Check sizes and get access pointers:
	assert(x.size()>0);
	assert(x[0].gInfo);
	const std::vector<double>& radialCoord = getRadialCoordinate(x[0]);
	std::vector<const double*> xData(x.size());
	for(unsigned i=0; i<x.size(); i++)
	{	assert(x[i].gInfo == x[0].gInfo);
		xData[i] = x[i].data();
	}
	//Write file:
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Could not open %s for writing.\n", filename)
	for(int j=0; j<x[0].gInfo->S; j++)
	{	fprintf(fp, "%le", radialCoord[j]);
		for(unsigned i=0; i<x.size(); i++)
			fprintf(fp, "\t%le", xData[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

#undef TCollection
#endif // FLUID1D_CORE1D_DATACOLLECTION_H
