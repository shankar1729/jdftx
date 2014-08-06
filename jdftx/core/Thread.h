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

#ifndef JDFTX_CORE_THREAD_H
#define JDFTX_CORE_THREAD_H

//! @addtogroup parallel
//! @{

/**
@file Thread.h
@brief Utilities for threading (wrappers around std::thread)
*/

#include <core/Util.h>
#include <thread>
#include <mutex>
#include <unistd.h>

extern int nProcsAvailable; //!< number of available processors (initialized to number of online processors, can be overriden)

/**
Operators should run multithreaded if this returns true,
and should run in a single thread if this returns false.
Note that all the thread launching functions in this file
automatically prevent nested threading, so operator codes
using those functions need not explicitly check this.

This only affects CPU threading, GPU operators should
only be called from a single thread anyway.
*/
bool shouldThreadOperators();

void suspendOperatorThreading(); //!< call from multi-threaded top-level code to disable threading within operators called from a parallel section
void resumeOperatorThreading(); //!< call after a parallel section in top-level code to resume threading within subsequent operator calls


/**
@brief A simple utility for running muliple threads

Given a callable object func and an argument list args, this routine calls nThreads threads
of func invoked as func(iMin, iMax, args). The nJobs jobs are evenly split between all the threads;
each instance of func should handle job index i satisfying iMin <= i < iMax.

If nJobs <= 0, the behaviour changes: the function is invoked as func(iThread, nThreads, args)
instead, where 0 <= iThread < nThreads. This mode allows for more flexible threading than the
evenly split job management indicated above. This could be used as a convenient interface for
launching threads for any parallel routine requiring as many threads as processors.

@param nThreads Number of threads to launch (if <=0, as many as processors on system)
@param func The function / object with operator() to invoke in a multithreaded fashion
@param nJobs The number of jobs to be split between the various func threads
@param args Arguments to pass to func
*/
template<typename Callable,typename ... Args>
void threadLaunch(int nThreads, Callable* func, size_t nJobs, Args... args);

/**
Same as threadLaunch(int nThreads, Callable* func, size_t nJobs, Args... args)
with nThreads = number of online processors.
*/
template<typename Callable,typename ... Args>
void threadLaunch(Callable* func, size_t nJobs, Args... args);


//! Maintain thread timing statistics and automatically choose the optimum number of threads
class AutoThreadCount
{
public:
	//! @param minThreads minimum number of threads to try
	//! @param minStats minimum number of tries for each thread count
	//! @param name if given, print debug messages & stats with that header
	//! @param fpLog debug logging stream
	AutoThreadCount(int minThreads=1, int minStats=3, const char* name=0, FILE* fpLog=stdout);
	~AutoThreadCount();
private:
	template<typename Callable,typename ... Args>
	friend void threadLaunch(AutoThreadCount*, Callable*, size_t, Args... args);

	int getThreadCount();
	void submitTime(int, double);
	double* time; int* count; int minThreads, nMax, nOpt, minStats; bool settled;
	char* name; FILE* fpLog;
};

/**
Same as threadLaunch(int nThreads, Callable* func, size_t nJobs, Args... args)
but with nThreads automatically adjusted over multiple runs using an AutoThreadCount object
If the AutoThreadCount pointer is null, as many threads as processors will be launched
*/
template<typename Callable,typename ... Args>
void threadLaunch(AutoThreadCount*, Callable* func, size_t nJobs, Args... args);


/**
@brief A parallelized loop

Given a callable object func and an argument list args, this routine calls func(i, args)
for each i in [0:nIter-1], and accumulates the return values.

Note that the calls are made from multiple threads, so func must
be thread safe. (Hint: pass mutexes as a part of args if synchronization
is required).

As many threads as online processors are launched and the nIter iterations are evenly split
between all the threads. Threaded loops will become single threaded if suspendOperatorThreading().

@param func The function / object with operator() to be looped over
@param nIter The number of loop 'iterations'
@param args Arguments to pass to func
*/
template<typename Callable,typename ... Args>
void threadedLoop(Callable* func, size_t nIter, Args... args);


/**
@brief A parallelized loop with an accumulated return value

Same as #threadedLoop, but func() returns double, which is summed over and returned
@return Accumulated return value of all calls to func()
*/
template<typename Callable,typename ... Args>
double threadedAccumulate(Callable* func, size_t nIter, Args... args);

//! @}


//###################################################################################################
//####  Implementation  ####
//##########################
//! @cond

template<typename Callable,typename ... Args>
void threadLaunch(int nThreads, Callable* func, size_t nJobs, Args... args)
{	if(nThreads<=0) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	if(nThreads>1) suspendOperatorThreading(); //Prevent func and anything it calls from launching nested threads
	std::thread** tArr = new std::thread*[nThreads-1];
	for(int t=0; t<nThreads; t++)
	{	size_t i1 = (nJobs>0 ? (  t   * nJobs)/nThreads : t);
		size_t i2 = (nJobs>0 ? ((t+1) * nJobs)/nThreads : nThreads);
		if(t<nThreads-1) tArr[t] = new std::thread(func, i1, i2, args...);
		else (*func)(i1, i2, args...);
	}
	for(int t=0; t<nThreads-1; t++)
	{	tArr[t]->join();
		delete tArr[t];
	}
	delete[] tArr;
	if(nThreads>1) resumeOperatorThreading(); //End nested threading guard section
}

template<typename Callable,typename ... Args>
void threadLaunch(Callable* func, size_t nJobs, Args... args)
{	threadLaunch(0, func, nJobs, args...);
}


template<typename Callable,typename ... Args>
void threadLaunch(AutoThreadCount* atc, Callable* func, size_t nJobs, Args... args)
{	if(!atc)
		threadLaunch(func, nJobs, args...);
	else
	{	int nThreads = atc->getThreadCount();
		double runTime = clock_us();
		threadLaunch(nThreads, func, nJobs, args...);
		runTime = clock_us() - runTime;
		atc->submitTime(nThreads, runTime);
	}
}


template<typename Callable,typename ... Args>
void threadedLoop_sub(size_t iMin, size_t iMax, Callable* func, Args... args)
{	for(size_t i=iMin; i<iMax; i++) (*func)(i, args...);
}
template<typename Callable,typename ... Args>
void threadedLoop(Callable* func, size_t nIter, Args... args)
{	threadLaunch(threadedLoop_sub<Callable,Args...>, nIter, func, args...);
}

template<typename Callable,typename ... Args>
void threadedAccumulate_sub(size_t iMin, size_t iMax, Callable* func, double* accumTot, std::mutex* m, Args... args)
{	register double accum=0.0;
	for(size_t i=iMin; i<iMax; i++) accum += (*func)(i, args...);
	m->lock(); *accumTot += accum; m->unlock();
}
template<typename Callable,typename ... Args>
double threadedAccumulate(Callable* func, size_t nIter, Args... args)
{	double accumTot=0.0;
	std::mutex m;
	threadLaunch(threadedAccumulate_sub<Callable,Args...>, nIter, func, &accumTot, &m, args...);
	return accumTot;
}

//!@endcond
#endif // JDFTX_CORE_THREAD_H
