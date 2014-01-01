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

#ifndef JDFTX_CORE_UTIL_H
#define JDFTX_CORE_UTIL_H

/** @file Util.h
@brief Miscellaneous utilities
*/

#include <core/MPIUtil.h>
#include <map>
#include <array>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>
#include <execinfo.h>

//------------- Common Initialization -----------------

extern bool killFlag; //!< Flag set by signal handlers - all compute loops should quit cleanly when this is set
extern MPIUtil* mpiUtil;
extern bool mpiDebugLog; //!< If true, all processes output to seperate debug log files, otherwise only head process outputs (set before calling initSystem())
void printVersionBanner(); //!< Print package name, version, revision etc. to log
void initSystem(int argc, char** argv); //!< Init MPI (if not already done), print banner, set up threads (play nice with job schedulers), GPU and signal handlers
void finalizeSystem(bool successful=true); //!< Clean-up corresponding to initSystem(), final messages (depending on successful) and clean-up MPI

//----------------- Profiling --------------------------

//! @brief Time of day in microseconds
inline double clock_us()
{	timeval tv;
	gettimeofday(&tv,NULL);
	return ((tv.tv_sec & 0x1fffff) * 1e6) + tv.tv_usec;
}

/** @brief Time a code section and print the timing (with title).\n
Code must be self-contained as a scope, as it will be surrounded by { }, and should not define "runTime"
@param title Name to use when printing the timing info for the code block
@param fp Stream to output the timing info to
@param code The code block to time (must be a self-contained scope)
*/
#define TIME(title,fp,code) \
{	double runTime = clock_us(); \
	{ code } \
	runTime = clock_us() - runTime; \
	fprintf(fp, "%s took %.2le s.\n", title, runTime*1e-6); \
}

//! Quick drop-in profiler for any function. Usage:
//! * Create a static object of this class in the function
//! * Call start and stop before and after the section to be timed
//! * Timing statistics of the code block will be printed on exit
#ifdef ENABLE_PROFILING
class StopWatch
{
public:
	StopWatch(string name);
	void start();
	void stop();
	void print() const;
private:
	double tPrev, Ttot, TsqTot; int nT;
	string name;
};
#else //ENABLE_PROFILING
//Version which does no profiling for release versions
class StopWatch
{
public:
	StopWatch(string name) {}
	void start() {}
	void stop() {}
};
#endif //ENABLE_PROFILING



// -----------  Debugging ---------------
void printStack(); //!< Print a minimal stack trace (convenient for debugging)
void gdbStackTraceExit(int code); //!< Exit on error with a more in-depth stack trace
int assertStackTraceExit(const char* expr, const char* function, const char* file, long line); //!< stack trace on failed assertions
//! A custom assertion with stack trace (NOTE: enabled in release modes as well)
#define assert(expr) \
	(void)((expr) ? 0 : assertStackTraceExit(#expr, __func__, __FILE__, __LINE__))


// -----------  Logging ---------------
extern FILE* globalLog;
extern FILE* nullLog; //!< pointer to /dev/null
void logSuspend(); //!< temporarily disable all log output (until logResume())
void logResume(); //!< re-enable logging after a logSuspend() call
#define logPrintf(...) fprintf(globalLog, __VA_ARGS__)
#define logFlush() fflush(globalLog)
#define die(...) \
	{	fprintf(globalLog, __VA_ARGS__); \
		if(mpiUtil->isHead() && globalLog != stdout) \
			fprintf(stderr, __VA_ARGS__); \
		finalizeSystem(false); \
		exit(1); \
	}


//--------------- Citations --------------------
namespace Citations
{	//!Add a citation to a paper with a reason
	//!(The same paper could be cited multiple times for different reasons)
	void add(string reason, string paper);
	
	//!Print the list of citations (with reasons) to the specified stream
	void print(FILE* fp=globalLog);
}


//--------------- Miscellaneous ------------------


#include <sys/stat.h>
//! Get the size of a file
off_t fileSize(const char *filename);

//! For any x and y>0, compute z = x % y such that 0 <= z < y
inline uint16_t positiveRemainder(int16_t x, uint16_t y)
{	register int16_t xMody = x % y;
	if(xMody < 0) return uint16_t(y + xMody);
	else return xMody;
}

//! Check if an integer is suitable for Fast Fourier Transforms (small prime factors only)
inline bool fftSuitable(int N)
{	static std::array<int,4> primes = {{2,3,5,7}};
	int tmpN = N;
	for(int p: primes) while(tmpN % p == 0) tmpN /= p;
	//All suitable prime factors taken out, we should be left with 1 for a suitable N:
	return (tmpN==1);
}

//! A template to ease option parsing (maps enums <--> strings)
template<typename Enum>
class EnumStringMap
{
	std::map<string,Enum> stringToEnum;
	std::map<Enum,string> enumToString;
	void addEntry() {}
	template<typename...Args> void addEntry(Enum e, const string& s, Args...args)
	{	stringToEnum[s]=e; enumToString[e]=s;
		addEntry(args...);
	}
public:
	//Initialize using the convenient syntax:
	// EnumStringMap mapname(enum1, string1, enum2, string2, ... as many as needed!)
	template<typename...Args> EnumStringMap(Args...args) { addEntry(args...); }

	//Match a key string to an enum, return false if not found
	bool getEnum(const char* key, Enum& e) const
	{	typename std::map<string,Enum>::const_iterator i = stringToEnum.find(key);
		if(i==stringToEnum.end()) return false;
		else
		{	e = i->second;
			return true;
		}
	}

	//Return the string representation of the enum:
	const char* getString(Enum e) const
	{	typename std::map<Enum,string>::const_iterator i = enumToString.find(e);
		return i->second.c_str();
	}

	string optionList() const
	{	typename std::map<string,Enum>::const_iterator i=stringToEnum.begin();
		string ret = i->first; i++;
		for(; i!=stringToEnum.end(); i++) ret += ("|"+i->first);
		return ret;
	}
};

#endif //JDFTX_CORE_UTIL_H
