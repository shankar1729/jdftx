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

#ifndef JDFTX_CORE_SCALARFIELDIO_H
#define JDFTX_CORE_SCALARFIELDIO_H

//! @addtogroup griddata
//! @{

/** @file ScalarFieldIO.h
@brief I/O utilities for the data arrays
*/

#include <core/ScalarField.h>
#include <core/vector3.h>
#include <core/Util.h>

#define Tptr std::shared_ptr<T>

//! Save the data in raw binary format to stream
template<typename T> void saveRawBinary(const Tptr& X, FILE* fp)
{	int nWrote = fwrite(X->data(), sizeof(typename T::DataType), X->nElem, fp);
	if(nWrote < X->nElem) die("Write failed after %d of %d records.\n", nWrote, X->nElem)
}
//! Save the data in raw binary format to file
template<typename T> void saveRawBinary(const Tptr& X, const char* filename)
{	FILE* fp = fopen(filename, "wb");
	if(!fp) die("Could not open '%s' for writing.\n", filename)
	saveRawBinary(X, fp);
	fclose(fp);
}

//! Load the data in raw binary format from stream
template<typename T> void loadRawBinary(Tptr& X, FILE* fp)
{	int nRead = fread(X->data(), sizeof(typename T::DataType), X->nElem, fp);
	if(nRead < X->nElem) die("Read failed after %d of %d records.\n", nRead, X->nElem)
}
//! Load the data in raw binary format from file
template<typename T> void loadRawBinary(Tptr& X, const char* filename)
{	FILE* fp = fopen(filename, "rb");
	if(!fp) die("Could not open '%s' for reading.\n", filename)
	
	off_t fLen = fileSize(filename);
	off_t expectedLen = sizeof(typename T::DataType) * X->nElem;
	if(fLen != expectedLen)
	{	die("\nLength of '%s' was %ld instead of the expected %ld bytes.\n"
				"Hint: Are you really reading the correct file?\n\n",
				filename, (unsigned long)fLen, (unsigned long)expectedLen);
	}
	
	loadRawBinary(X, fp);
	fclose(fp);
}

#undef Tptr

/** Save data to a raw binary along with a DataExplorer header
@param filenamePrefix Binary data is saved to filenamePrefix.bin with DataExplorer header filenamePrefix.dx
*/
void saveDX(const ScalarField&, const char* filenamePrefix);

/** Spherically average scalar fields about an arbitrary center (with Wigner-Seitz wrapping)
@param dataR The data to sphericalize and save
@param nColumns Number of ScalarField's in dataR[]
@param drFac is the spacing in radius as a fraction of the diameter of the sample box (R ./ S) (drFac << 1 is likely to give noisy results, particularly close to r=0)
@param center The origin for spherical coordinates [default = center of box (if null pointer is passed)]
@return The first array contains the radial grid, and the subsequent ones the spherically-averaged results, one for each dataR, and the last column contains the weight of the radial grid point
*/
std::vector< std::vector<double> > sphericalize(const ScalarField* dataR, int nColumns, double drFac=1.0, vector3<>* center=0);

/** Saves an array of real space data pointers to a multicolumn 1D 'sphericalized' file (for gnuplot)
@param dataR The data to sphericalize and save
@param nColumns Number of ScalarField's in dataR[]
@param filename Output file in which column 1 will be the radius, column 2 to nColumns+1 would be the sphericalized versions of dataR[0 to nColumns-1]
@param drFac is the spacing in radius as a fraction of the diameter of the sample box (R ./ S) (drFac << 1 is likely to give noisy results, particularly close to r=0)
@param center The origin for spherical coordinates [default = center of box (if null pointer is passed)]
*/
void saveSphericalized(const ScalarField* dataR, int nColumns, const char* filename, double drFac=1.0, vector3<>* center=0);


/** Saves an array of reciprocal space data pointers to a multicolumn 1D 'sphericalized' file (for gnuplot)
@param dataG The data to sphericalize (about G=0) and save
@param nColumns Number of ScalarFieldTilde's in dataG[]
@param filename Output file in which column 1 will be the radius, column 2 to nColumns+1 would be the sphericalized versions of dataG[0 to nColumns-1]
@param dGFac is the spacing in radius as a fraction of the diameter of the Brillouin zone (dGFac << 1 is likely to give noisy results, particularly close to G=0)
*/
void saveSphericalized(const ScalarFieldTilde* dataG, int nColumns, const char* filename, double dGFac=1.0);

//! @}
#endif // JDFTX_CORE_SCALARFIELDIO_H
