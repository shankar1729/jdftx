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

#ifndef JDFTX_CORE_DATACOLLECTION_H
#define JDFTX_CORE_DATACOLLECTION_H


//! @file DataCollection.h
//! @brief classes DataRptrCollection, DataGptrCollection and just enough operators to enable CG w.r.t to them

#include <core/Operators.h>
#include <vector>

#include <cstdio>

typedef std::vector<DataRptr> DataRptrCollection; //!< dynamic size collection of real space scalar fields
typedef std::vector<DataGptr> DataGptrCollection; //!< dynamic size collection of reciprocal space scalar fields
#define TptrCollection std::vector<std::shared_ptr<T> > //!< shorthand for templates below (undef'd at end of file)

//! Scale
template<typename T> TptrCollection& operator*=(TptrCollection& x, double alpha)
{ 	for(unsigned i=0; i<x.size(); i++) x[i] *= alpha;
	return x;
}

//! y += alpha x
template<typename T> void axpy(double alpha, const TptrCollection& x, TptrCollection& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) axpy(alpha, x[i], y[i]);
}

//! Inner product
template<typename T> double dot(const TptrCollection& x, const TptrCollection& y)
{	assert(x.size()==y.size());
	double ret = 0.0;
	for(unsigned i=0; i<x.size(); i++) ret += dot(x[i], y[i]);
	return ret;
}

//! Create a copy of the data (note operator= references same data since Tptr's are pointers!)
template<typename T> TptrCollection clone(const TptrCollection& x)
{	TptrCollection ret(x.size());
	for(unsigned i=0; i<x.size(); i++) ret[i] = clone(x[i]);
	return ret;
}

//! Initialize (non-null) data to zero
template<typename T> void initZero(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initZero(x[i]);
}

//! Allocate all null components and initialize them to zero
//! If N is non-zero, resize collection to N scalar fields
template<typename T> void nullToZero(TptrCollection& x, const GridInfo& gInfo, int N=0)
{	if(N) x.resize(N);
	for(unsigned i=0; i<x.size(); i++) nullToZero(x[i], gInfo);
}

//! Initialize to random numbers (uniform on 0 to 1)
template<typename T> void initRandomFlat(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandomFlat(x[i]);
}
//! Initialize to normal random numbers:
template<typename T> void randomize(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandom(x[i], 3.0);
}

template<typename T> void loadFromFile(TptrCollection& x, const char* filename)
{	//Checks for the correct filesize
	off_t expectedLen = 0;
	for(unsigned i=0; i<x.size(); i++){expectedLen += sizeof(typename T::DataType) * x[i]->nElem;}
	off_t fLen = fileSize(filename);
	if(fLen != expectedLen)
	{	die("\nLength of '%s' was %ld instead of the expected %ld bytes.\n"
				"Hint: Are you really reading the correct file?\n\n",
				filename, fLen, expectedLen);
	}
	
	FILE* fp = fopen(filename, "rb");
	if(!fp) die("Could not open %s for reading.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in loadFromFile(x,\"%s\").\n", i, filename)
		if(fread(x[i]->data(), sizeof(typename T::DataType), x[i]->nElem, fp) < unsigned(x[i]->nElem))
			die("File ended too soon while reading x[%d] in loadFromFile(x,\"%s\").\n", i, filename)
	}
	fclose(fp);
}

template<typename T> void saveToFile(const TptrCollection& x, const char* filename)
{	FILE* fp = fopen(filename, "wb");
	if(!fp) die("Could not open %s for writing.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in saveToFile(x,\"%s\").\n", i, filename)
		fwrite(x[i]->data(), sizeof(typename T::DataType), x[i]->nElem, fp);
	}
	fclose(fp);
}

#undef TptrCollection
#endif // JDFTX_CORE_DATACOLLECTION_H
