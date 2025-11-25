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
#include <functional>
#include <core/string.h>
#include <core/vector3.h>

//! Data structure to receive N-D buffers from Python (e.g. numpy.ndarray or torch.Tensor)
//! Keeping this double-precision real only for now; can template this later if needed.
struct NDarray
{	double* data;
	std::vector<size_t> shape;
	std::vector<size_t> strides;
	
	inline double& operator[](const std::vector<size_t>& index)
	{	size_t offset = 0;
		for(size_t dim=0; dim<shape.size(); dim++)
			offset += index[dim] * strides[dim];
		return data[offset];
	}
};

typedef std::function<void(NDarray, std::vector<NDarray>)> setCavity_t;


//! Wrapper that maximally isolates JDFTx internals from the python interface
class JDFTxWrapper
{
public:
	JDFTxWrapper(std::vector<std::pair<string, string>> inputs, bool variableCell, setCavity_t setCavity[2]);
	void minimize();
	void move(NDarray delta_positions, NDarray delta_R);
	void dumpEnd() const;
	double getEnergy() const;
	NDarray getForces() const;
	NDarray getStress() const;
	static std::vector<string> getCommands();
	static string getCommandDoc(string name);
private:
	bool variableCell;
	void compute();
	std::shared_ptr<class Everything> e;
	std::shared_ptr<class IonicMinimizer> imin;
	std::shared_ptr<class LatticeMinimizer> lmin;
	std::vector<vector3<>> forces; //flat array of Cartesian forces (converted from SpeciesInfo::forces)
};

#endif //JDFTX_ELECTRONIC_JDFTX_WRAPPER_H
