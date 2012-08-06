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

#ifndef JDFTX_COMMAND_PARAMLIST_H
#define JDFTX_COMMAND_PARAMLIST_H

#include <core/Util.h>

//! List of parameters for the command
class ParamList
{
	istringstream iss;
public:
	//! Construct given the string with all the parameters
	ParamList(string params) : iss(params) {}

	//! Retreive a paremeter of type T in t, with default value tDefault if the stream ends
	//! paramName is used to report an error if conversion failed
	//! If required=true, error will be reported instead of setting the default value on eof()
	template<typename T>
	void get(T& t, T tDefault, string paramName, bool required=false)
	{	iss.clear(); //clear previous errors
		iss >> t;
		if(iss.bad()) throw string("I/O error while reading parameter <"+paramName+">.");
		if(iss.eof())
		{	if(required) throw string("Parameter <"+paramName+"> must be specified.");
			else { t = tDefault; return; }
		}
		if(iss.fail()) throw string("Conversion of parameter <"+paramName+"> failed.");
	}

	//! Retreive an enum paremeter of type T in t, with default value tDefault if the stream ends
	//! tMap is used to convert the string to the enum
	//! paramName is used to report an error if conversion failed
	//! If required=true, error will be reported instead of setting the default value on eof()
	template<typename T>
	void get(T& t, T tDefault, EnumStringMap<T>& tMap, string paramName, bool required=false)
	{	iss.clear(); //clear previous errors
		string key;
		iss >> key;
		if(iss.bad()) throw string("I/O error while reading parameter <"+paramName+">.");
		if(iss.eof())
		{	if(required) throw string("Parameter <"+paramName+"> must be specified.");
			else { t = tDefault; return; }
		}
		if(!tMap.getEnum(key.c_str(), t))
			throw string("Parameter <"+paramName+"> must be one of "+tMap.optionList());
	}
};

#endif //JDFTX_COMMAND_PARAMLIST_H
