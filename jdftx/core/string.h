/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_ISTRING_H
#define JDFTX_CORE_ISTRING_H

//! @file string.h STL strings and streams with case insensitive comparison

#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>

//! Case insensitive character trait
struct ichar_traits : public std::char_traits<char>
{	static bool eq( char c1, char c2 ) { return toupper(c1) == toupper(c2); }
	static bool ne( char c1, char c2 ) { return toupper(c1) != toupper(c2); }
	static bool lt( char c1, char c2 ) { return toupper(c1) <  toupper(c2); }
	static int compare( const char* s1, const char* s2, size_t n ) { return strncasecmp( s1, s2, n ); }
	static const char* find( const char* s, int n, char a )
	{	for(int i=0; i<n; i++) if(toupper(s[i])==toupper(a)) return s+i;
		return 0;
	}
};
typedef std::basic_string<char, ichar_traits> string; //!< Case-insensitive string

//! Remove leading and trailing spaces from a string
inline void trim(string& s)
{	//Trim trailing whitespace:
	size_t endNonWS = s.find_last_not_of(" \t\n\r");
	if(endNonWS != string::npos)
		s.erase(endNonWS+1);
	//Trim leading whitespace:
	s.erase(0, s.find_first_not_of(" \t\n\r")); //deletes entire line if all whitespace
}

using std::istream;
using std::ostream;

inline istream& operator>>(istream& is, string& str) { std::string tmp; is >> tmp; str.assign(tmp.c_str()); return is; }
inline ostream& operator<<(ostream& os, const string& str) { os << str.c_str(); return os; }
inline istream& getline (istream& is, string& str, char delim='\n')
{	std::string tmp;
	getline(is, tmp, delim);
	str.assign(tmp.c_str());
	return is; 
}

struct ifstream : public std::ifstream
{	ifstream() {}
	explicit ifstream(const string& fname) : std::ifstream(fname.c_str()) {}
	void open(const string& fname) { std::ifstream::open(fname.c_str()); }
};
struct ofstream : public std::ofstream
{	ofstream() {}
	explicit ofstream(const string& fname) : std::ofstream(fname.c_str()) {}
	void open(const string& fname) { std::ofstream::open(fname.c_str()); }
};

struct istringstream : public std::istringstream
{	istringstream() {}
	explicit istringstream(const string& s) { std::istringstream::str(std::string(s.c_str())); }
	void str(const string& s) { std::istringstream::str(std::string(s.c_str())); }
};
struct ostringstream : public std::ostringstream
{	string str() const { return string(std::ostringstream::str().c_str()); }
};

#endif
