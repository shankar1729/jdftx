/*-------------------------------------------------------------------
Copyright 2025 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_JDFTX_WRAPPER_H
#define JDFTX_ELECTRONIC_JDFTX_WRAPPER_H

#include <memory>
#include <vector>
#include <core/string.h>

//! Data structure to receive N-D buffers from Python (e.g. numpy.ndarray or torch.Tensor)
//! Keeping this double-precision real only for now; can template this later if needed.
struct NDarray
{	double* data;
	std::vector<size_t> shape;
	std::vector<size_t> strides;
	
	inline size_t getOffset(const std::vector<size_t>& index) const
	{	size_t offset = 0;
		for(size_t dim=0; dim<shape.size(); dim++)
			offset += index[dim] * strides[dim];
		return offset;
	}
	inline double& operator[](const std::vector<size_t>& index) {return data[getOffset(index)];}
	inline const double& operator[](const std::vector<size_t>& index) const {return data[getOffset(index)];}
};


//! Wrapper that maximally isolates JDFTx internals from the python interface
class JDFTxWrapper
{
public:
	JDFTxWrapper(std::vector<std::pair<string, string>> inputs, bool variableCell);
	void minimize();
	void move(NDarray delta_positions, NDarray delta_R);
private:
	bool variableCell;
	std::shared_ptr<class Everything> e;
	std::shared_ptr<class IonicMinimizer> imin;
	std::shared_ptr<class LatticeMinimizer> lmin;
};

#endif //JDFTX_ELECTRONIC_JDFTX_WRAPPER_H
