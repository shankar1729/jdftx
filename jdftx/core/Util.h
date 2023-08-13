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

//! @addtogroup Utilities
//! @{

//! @file Util.h Miscellaneous utilities

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

extern string inputBasename; //!< Basename of input file or "stdin" that can be used as a default run-name
extern bool killFlag; //!< Flag set by signal handlers - all compute loops should quit cleanly when this is set
extern MPIUtil* mpiWorld; //!< MPI across all processes
extern MPIUtil* mpiGroup; //!< MPI within current group of processes
extern MPIUtil* mpiGroupHead; //!< MPI across equal ranks in each group
extern bool mpiDebugLog; //!< If true, all processes output to seperate debug log files, otherwise only head process outputs (set before calling initSystem())
extern size_t mempoolSize; //!< If non-zero, size of memory pool managed internally by JDFTx

//! Parameters used for common initialization functions
struct InitParams
{	//Input parameters:
	const char* description; //!< description of program used when printing usage
	class Everything* e; //!< pointer to use when calling template
	InitParams(const char* description=0, class Everything* e=0);
	//Output parameters retrieved from command-line:
	string inputFilename; //!< name of input file
	bool dryRun; //!< whether this is a dry run
	bool printDefaults; //!< whether to print default commands
	//Optional parameters useful when calling from outside JDFTx:
	const char* packageName; //!< package name dispalyed in banner
	const char* versionString; //!< version string displayed in banner
	const char* versionHash; //!< git hash displayed in banner (if any)
};
void printVersionBanner(const InitParams* ip=0); //!< Print package name, version, revision etc. to log
void initSystem(int argc, char** argv, const InitParams* ip=0); //!< Init MPI (if not already done), print banner, set up threads (play nice with job schedulers), GPU and signal handlers
void initSystemCmdline(int argc, char** argv, InitParams& ip); //!< initSystem along with commandline options
void finalizeSystem(bool successful=true); //!< Clean-up corresponding to initSystem(), final messages (depending on successful) and clean-up MPI

//----------------- Profiling --------------------------

double clock_us(); //! @brief Elapsed time in microseconds (from start of program)
double clock_sec(); //! @brief Elapsed time in seconds (from start of program)

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
void printStack(bool detailedStackScript=false); //!< Print a minimal stack trace and optionally write a script that, when run, will print a more detailed stacktrace
void stackTraceExit(int code); //!< Exit on error with stack trace
int assertStackTraceExit(const char* expr, const char* function, const char* file, long line); //!< stack trace on failed assertions
//! A custom assertion with stack trace (NOTE: enabled in release modes as well)
#define assert(expr) \
	(void)((expr) ? 0 : assertStackTraceExit(#expr, __func__, __FILE__, __LINE__))


// -----------  Logging ---------------
extern FILE* globalLog;
extern FILE* nullLog; //!< pointer to /dev/null
void logSuspend(); //!< temporarily disable all log output (until logResume())
void logResume(); //!< re-enable logging after a logSuspend() call

#define logPrintf(...) fprintf(globalLog, __VA_ARGS__) //!< printf() for log files
#define logFlush() fflush(globalLog) //!< fflush() for log files

//! @brief Quit with an error message (formatted using printf()). Must be called from all processes.
#define die(...) \
	{	fprintf(globalLog, __VA_ARGS__); \
		if(mpiWorld->isHead() && globalLog != stdout) \
			fprintf(stderr, __VA_ARGS__); \
		finalizeSystem(false); \
		exit(1); \
	}

//! @brief Version of die that should only be used when it is impossible to guarantee synchronized calls from all processes
#define die_alone(...) \
	{	fprintf(globalLog, __VA_ARGS__); \
		fflush(globalLog); \
		if(mpiWorld->isHead() && globalLog != stdout) \
			fprintf(stderr, __VA_ARGS__); \
		mpiWorld->exit(1); \
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

#include <inttypes.h>
#ifndef PRIdPTR
	#define PRIdPTR "zd" //For pre-C++11 compilers
#endif

//Endianness utilities (all binary I/O is from little-endian files regardless of operating endianness):
void convertToLE(void* ptr, size_t size, size_t nmemb); //!< Convert data from operating endianness to little-endian
void convertFromLE(void* ptr, size_t size, size_t nmemb); //!< Convert data from little-endian to operating endianness
size_t freadLE(void *ptr, size_t size, size_t nmemb, FILE* fp); //!< Read from a little-endian binary file, regardless of operating endianness
size_t fwriteLE(const void *ptr, size_t size, size_t nmemb, FILE *fp); //!< Write to a little-endian binary file, regardless of operating endianness

//! For any x and y>0, compute z = x % y such that 0 <= z < y
inline uint16_t positiveRemainder(int16_t x, uint16_t y)
{	int16_t xMody = x % y;
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

//! @}
#endif //JDFTX_CORE_UTIL_H
