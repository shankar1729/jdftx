/*-------------------------------------------------------------------
Copyright 2017 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_H5IO_H
#define JDFTX_CORE_H5IO_H

//! @addtogroup Output
//! @{

/** @file H5io.h
@brief HDF5 helper routines
*/

#include <core/Util.h>
#include <array>

#ifdef HDF5_ENABLED
#include "hdf5.h"

inline hid_t h5createGroup(hid_t parent, const char* name);
template<typename T> void h5writeScalar(hid_t fid, const char* dname, const T& data); //Write scalar to a rank-0 dataset
template<typename T> void h5writeVector(hid_t fid, const char* dname, const std::vector<T>& data); //Collectively write contiguous array to a 1D dataset when all the data is available on all the processes
template<typename T> void h5writeVector(hid_t fid, const char* dname, const T* data, hsize_t nData); //Collectively write contiguous array to a 1D dataset when all the data is available on all the processes
template<typename T> void h5writeVector(hid_t fid, const char* dname, const T* data, const hsize_t* dims, hsize_t rank); //Collectively write contiguous array to a nD dataset when all the data is available on all the processes

//! @}

//------------ Implementations -------------------
//!@cond

template<typename T> struct h5type;
template<> struct h5type<int> { static hid_t get() { return H5T_NATIVE_INT; } };
template<> struct h5type<double> { static hid_t get() { return H5T_NATIVE_DOUBLE; } };

inline hid_t h5createGroup(hid_t parent, const char* name)
{	hid_t gid = H5Gcreate(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(gid<0) die("Error creating group '%s' in HDF5 file.\n", name);
	return gid;
}

template<typename T> void h5writeScalar(hid_t fid, const char* dname, const T& data)
{	hid_t dataType = h5type<T>::get();
	//Create dataset:
	hid_t sid = H5Screate(H5S_SCALAR);
	hid_t did = H5Dcreate(fid, dname, dataType, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sid);
	if(did<0) die("Could not create dataset '%s' in HDF5 file.\n", dname);
	//Write data:
	H5Dwrite(did, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
	H5Dclose(did);
}

template<typename T> void h5writeVector(hid_t fid, const char* dname, const std::vector<T>& data)
{	h5writeVector(fid, dname, data.data(), data.size());
}
template<typename T> void h5writeVector(hid_t fid, const char* dname, const T* data, hsize_t nData)
{	hsize_t dims[1] = { nData };
	h5writeVector(fid, dname, data, dims, 1);
}
template<typename T> void h5writeVector(hid_t fid, const char* dname, const T* data, const hsize_t* dims, hsize_t rank)
{	hid_t dataType = h5type<T>::get();
	//Create dataset:
	hid_t sid = H5Screate_simple(rank, dims, NULL);
	hid_t did = H5Dcreate(fid, dname, dataType, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sid);
	if(did<0) die("Could not create dataset '%s' in HDF5 file.\n", dname);
	//Create collective write property:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Write data:
	H5Dwrite(did, dataType, H5S_ALL, H5S_ALL, plid, data);
	//Cleanup:
	H5Dclose(did);
	H5Pclose(plid);
}

//!@endcond
#endif //HDF5_ENABLED
#endif //JDFTX_CORE_H5IO_H
