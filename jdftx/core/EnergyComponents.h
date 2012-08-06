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

#ifndef JDFTX_CORE_ENERGYCOMPONENTS_H
#define JDFTX_CORE_ENERGYCOMPONENTS_H

//! @file EnergyComponents.h
//! @brief Represent components of the (free) energy

#include <map>
#include <string>
#include <cstdio>

/** Simply a map of named components with a few extra features!
	Proper use: E["name"] = value; or E["name"] += value
	Casting to a double returns the total.
	Assignment or increment by a complex happens in a "default" component name
*/
class EnergyComponents : public std::map<string,double>
{
public:

	//! Set to a simple complex
	void operator=(const double& value)
	{	clear(); //remove all named components
		(*this)["default"] = value; //set the value to a default channel
	}

	//! Accumulate a simple complex
	void operator+=(const double& value)
	{	(*this)["default"] += value; //accumulate to the default channel
	}

	//! Return sum of all components when cast to a simple number
	// Allows using templated CG without modification!
	operator double() const
	{	double ret=0.0;
		for(const_iterator i=begin(); i!=end(); i++)
		{	ret += i->second;
		}
		return ret;
	}

	//! Print to a stream:
	//! @param fp Output stream
	//! @param nonzeroOnly If true, print only non-zero components
	//! @param format Output format for each component - must have one string and one floating point conversion spec each
	void print(FILE* fp, bool nonzeroOnly, const char* format="\t%s = %le\n")
	{	for(const_iterator i=begin(); i!=end(); i++)
			if(!nonzeroOnly || i->second!=0.0)
				fprintf(fp, format, i->first.c_str(), i->second);
	}
};

#endif //JDFTX_CORE_ENERGYCOMPONENTS_H
