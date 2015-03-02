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

/** @file ParamList.h
@brief Helper class for parsing command lines in input file
*/

//! @brief Wrapper to std::istringstream that eases parsing of input file command lines
class ParamList
{
	istringstream iss;
public:
	//! Construct given the string with all the parameters
	ParamList(string params) : iss(params) {}

	//! Rewind to beginning of stream (useful for commands with multiple alternate syntaxes)
	void rewind() { iss.seekg(0, std::ios_base::beg); }
	
	//! @brief Retrieve a parameter from stream or, if stream ends, set it to a default value.
	//! Note that t is set to tDefault if any error occurs.
	//! @tparam T data type of parameter, which could be anything that std::istream supports
	//! @param t Destination for the retrieved parameter
	//! @param tDefault Default vaue of the parameter to set if stream ends while parsing
	//! @param paramName Name of the parameter to use in error messages
	//! @param required If true, report an error if the stream ends instead of using the default.
	template<typename T>
	void get(T& t, T tDefault, string paramName, bool required=false)
	{	iss.clear(); //clear previous errors
		iss >> t;
		if(iss.bad()) throw string("I/O error while reading parameter <"+paramName+">.");
		if(iss.eof())
		{	t = tDefault;
			if(required) throw string("Parameter <"+paramName+"> must be specified.");
			else return;
		}
		if(iss.fail()) { t = tDefault; throw string("Conversion of parameter <"+paramName+"> failed."); }
	}

	//! @brief Retreive an enum parameter from stream with default handling 
	//! @tparam T Type of the enum parameter
	//! @param t Destination for the retrieved parameter
	//! @param tDefault Default vaue of the parameter to set if stream ends while parsing
	//! @param tMap Map between the enum values and corresponding string names
	//! @param paramName Name of the parameter to use in error messages
	//! @param required If true, report an error if the stream ends instead of using the default.
	template<typename T>
	void get(T& t, T tDefault, const EnumStringMap<T>& tMap, string paramName, bool required=false)
	{	iss.clear(); //clear previous errors
		string key;
		iss >> key;
		if(iss.bad()) throw string("I/O error while reading parameter <"+paramName+">.");
		if(iss.eof())
		{	t = tDefault;
			if(required) throw string("Parameter <"+paramName+"> must be specified.");
			else return;
		}
		if(!tMap.getEnum(key.c_str(), t))
		{	t = tDefault;
			throw string("Parameter <"+paramName+"> must be one of "+tMap.optionList());
		}
	}
	
	//!Get the section of the input string not yet parsed
	string getRemainder()
	{	if(iss.eof()) return string();
		size_t curPos = iss.tellg();
		iss.seekg(0, std::ios_base::end);
		size_t endPos = iss.tellg();
		if(endPos > curPos)
		{	string buf(endPos-curPos, 0);
			iss.seekg(curPos);
			iss.read(&buf.at(0), buf.length());
			trim(buf);
			return buf;
		}
		else return string();
	}
};

#endif //JDFTX_COMMAND_PARAMLIST_H
