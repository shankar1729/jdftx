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

#include <core/Util.h>
#include <core/Thread.h>
#include <cmath>
#include <csignal>
#include <list>
#include <vector>

#ifndef __APPLE__
#include <sys/prctl.h>
#endif

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#endif

#include <config.h> //This file is generated during build based on SVN version etc.

// Get the size of a file
off_t fileSize(const char *filename)
{	struct stat st;
	if(stat(filename, &st) == 0) return st.st_size;
    return -1;
}

//Program banner
void printVersionBanner()
{	logPrintf("\n*************** " PACKAGE_NAME " " VERSION_MAJOR_MINOR_PATCH
		" (svn revision " VERSION_SVN_REVISION ") " PACKAGE_DESCRIPTION
		" ****************\n\n");
}

//Signal handlers:
bool killFlag = false; //!< The signal handlers set this flag to suggest clean termination
void sigIntHandler(int sig); //!< Handle Ctrl+C (quit cleanly, abruptly or not at all - based on user input)
void sigQuitHandler(int sig); //!< Handle Ctrl+\ (cleanly quit after current iteration, without prompting)
void sigErrorHandler(int sig); //!< Exit with a stack trace on segfaults and aborts (to ease debugging)
void resetHandlers()
{	signal(SIGINT, SIG_IGN);
	signal(SIGQUIT, SIG_IGN);
}
void registerHandlers()
{	signal(SIGINT, sigIntHandler);
	signal(SIGQUIT, sigQuitHandler);
	signal(SIGSEGV, sigErrorHandler);
	signal(SIGABRT, sigErrorHandler);
}
void sigIntHandler(int sig)
{	if(feof(stdin)) mpiUtil->exit(0);
	resetHandlers();
	printf(
		"\n---------------------------------------------\n"
		"Received SIGINT (Ctrl+C), do you want to:\n"
		"\t[Q] Quit rightaway?\n"
		"\t[A] Quit cleanly after current iteration?\n"
		"\t[I] Ignore and continue normally?\n");
	while(1)
	{
		printf("Enter [Q/A/I]: "); fflush(stdout);
		char c = getchar();
		switch(c)
		{	case 'q': case 'Q': printf("Quitting now ...\n"); mpiUtil->exit(0);
			case 'a': case 'A':
				printf("Will quit after current iteration ...\n");
				killFlag = true; registerHandlers(); return;
			case 'i': case 'I':
				printf("Ignoring and continuing normally ...\n");
				registerHandlers(); return;
			default:
				printf("Unrecognized input.\n");
		}
	}
}
void sigQuitHandler(int sig)
{	resetHandlers();
	logPrintf("\n"
		"+---------------------------------------------------------------+\n"
		"|  Received SIGQUIT, will quit cleanly after current iteration  |\n"
		"+---------------------------------------------------------------+\n\n");
	logFlush();
	killFlag = true;
	registerHandlers();
}
void sigErrorHandler(int sig)
{	fprintf(stderr, sig==SIGSEGV ? "Segmentation Fault.\n" : "Aborted.\n");
	gdbStackTraceExit(1);
}

FILE* globalLog = stdout; // this might be replaced by a file pointer in main before calling initSystem()
FILE* globalLogOrig;
FILE* nullLog = 0;

void logSuspend()
{	if(nullLog) globalLog = nullLog;
}

void logResume()
{	globalLog = globalLogOrig;
}

MPIUtil* mpiUtil = 0;
bool mpiDebugLog = false;
static double startTime_us; //Time at which system was initialized in microseconds

void initSystem(int argc, char** argv)
{
	if(!mpiUtil) mpiUtil = new MPIUtil(argc, argv);
	nullLog = fopen("/dev/null", "w");
	if(!mpiUtil->isHead())
	{	if(mpiDebugLog)
		{	char fname[256]; sprintf(fname, "jdftx.%d.mpiDebugLog", mpiUtil->iProcess());
			globalLog = fopen(fname, "w");
		}
		else globalLog = nullLog;
	}
	globalLogOrig = globalLog;
	
	//Print a welcome banner with useful information
	printVersionBanner();
	time_t startTime = time(0);
	startTime_us = clock_us();
	logPrintf("Start date and time: %s", ctime(&startTime)); //note ctime output has a "\n" at the end
	char hostname[256]; gethostname(hostname, 256);
	logPrintf("Running on host: %s\n", hostname);
	logPrintf("Executable %s with ", argv[0]);
	if(argc>1)
	{	logPrintf("command-line:");
		for(int i=1; i<argc; i++) logPrintf(" %s", argv[i]);
		logPrintf("\n");
	}
	else logPrintf(" empty command-line.\n");
	
	registerHandlers();
	
	#ifdef GPU_ENABLED
	if(!gpuInit(globalLog)) die("gpuInit() failed\n\n")
	#endif
	
	//Limit thread count if running within SLURM:
	const char* slurmCpusPerTask = getenv("SLURM_CPUS_PER_TASK");
	if(slurmCpusPerTask)
	{	int nThreadsMax;
		if(sscanf(slurmCpusPerTask, "%d", &nThreadsMax)==1)
			nProcsAvailable = nThreadsMax; //Slurm spec found, update available processor count (Thread.h)
		else
			logPrintf("Could not determine thread count from SLURM_CPUS_PER_TASK=\"%s\".\n", slurmCpusPerTask);
	}
	logPrintf("Current process will run with a maximum of %d cpu threads.\n", nProcsAvailable);
	
	//Print total resources used by run:
	{	int nResources[2] = { nProcsAvailable, 0 };
		#ifdef GPU_ENABLED
		nResources[1] = 1;
		#endif
		mpiUtil->allReduce(nResources, 2, MPIUtil::ReduceSum);
		logPrintf("Run totals: %d processes, %d threads, %d GPUs\n", mpiUtil->nProcesses(), nResources[0], nResources[1]);
	}
	
	//Add citations to the code and general framework:
	Citations::add("Software package", "R. Sundararaman, K. Letchworth-Weaver and T.A. Arias, JDFTx, available from http://jdftx.sourceforge.net (2012)");
	Citations::add("Algebraic framework", "S. Ismail-Beigi and T.A. Arias, Computer Physics Communications 128, 1 (2000)");
}


#ifdef ENABLE_PROFILING
void stopWatchManager(const StopWatch* addWatch=0)
{	static std::vector<const StopWatch*> watches; //static array of all watches
	if(addWatch) watches.push_back(addWatch);
	else //print timings:
	{	logPrintf("\n");
		for(const StopWatch* watch: watches) watch->print();
	}
}
#endif // ENABLE_PROFILING


void finalizeSystem(bool successful)
{
	time_t endTime = time(0);
	char* endTimeString = ctime(&endTime);
	endTimeString[strlen(endTimeString)-1] = 0; //get rid of the newline in output of ctime
	double durationSec = 1e-6*(clock_us() - startTime_us);
	int durationDays = floor(durationSec/86400.); durationSec -= 86400.*durationDays;
	int durationHrs  = floor(durationSec/3600.);  durationSec -= 3600.*durationHrs;
	int durationMin  = floor(durationSec/60.);    durationSec -= 60.*durationMin;
	logPrintf("End date and time: %s  (Duration: %d-%d:%02d:%05.2lf)\n",
		endTimeString, durationDays, durationHrs, durationMin, durationSec);
	
	if(successful) logPrintf("Done!\n");
	else
	{	logPrintf("Failed.\n");
		if(mpiUtil->isHead() && globalLog != stdout)
			fprintf(stderr, "Failed.\n");
	}
	
	#ifdef ENABLE_PROFILING
	stopWatchManager();
	#endif
	
	if(!mpiUtil->isHead())
	{	if(mpiDebugLog) fclose(globalLog);
		globalLog = 0;
	}
	fclose(nullLog);
	delete mpiUtil;
}


#ifdef ENABLE_PROFILING
StopWatch::StopWatch(string name) : Ttot(0), TsqTot(0), nT(0), name(name) { stopWatchManager(this); }
void StopWatch::start()
{
	#ifdef GPU_ENABLED
	cudaThreadSynchronize();
	#endif
	tPrev = clock_us();
}
void StopWatch::stop()
{
	#ifdef GPU_ENABLED
	cudaThreadSynchronize();
	#endif
	double T = clock_us()-tPrev;
	Ttot+=T; TsqTot+=T*T; nT++;
}
void StopWatch::print() const
{	if(nT)
	{	double meanT = Ttot/nT;
		double sigmaT = sqrt(TsqTot/nT - meanT*meanT);
		logPrintf("PROFILER: %30s %12.6lf +/- %12.6lf s, %4d calls, %13.6lf s total\n",
			name.c_str(), meanT*1e-6, sigmaT*1e-6, nT, Ttot*1e-6);
	}
}
#endif //ENABLE_PROFILING


// Print a minimal stack trace (convenient for debugging)
void printStack()
{	void* tracePtrs[100];
	int count = backtrace(tracePtrs, 100);
	char** funcNames = backtrace_symbols(tracePtrs, count);
	for(int i=0; i<count; i++) printf("\t%s\n", funcNames[i]);
	free(funcNames);
}

// Exit on error with a more in-depth stack trace
void gdbStackTraceExit(int code)
{	// From http://stackoverflow.com/questions/4636456/stack-trace-for-c-using-gcc/4732119#4732119
	char pid_buf[30]; sprintf(pid_buf, "%d", getpid());
    char name_buf[512]; name_buf[readlink("/proc/self/exe", name_buf, 511)]=0;
	int fdPipe[2]; if(pipe(fdPipe)) { printf("Error creating pipe.\n"); mpiUtil->exit(code); }
	char message = '\n'; //some random character for sync
    int child_pid = fork();
    if(!child_pid)
	{	dup2(2,1); //redirect output to stderr
		//Wait for ptrace permissions to be set by parent:
		close(fdPipe[1]);
		while(!read(fdPipe[0], &message, 1));
		close(fdPipe[0]);
		//Attach gdb:
		fprintf(stdout,"\n\nStack trace:\n");
		execlp("gdb", "gdb", "--batch", "-n", "-ex", "bt", name_buf, pid_buf, NULL);
		abort(); //If gdb failed to start
    }
    else
	{	
		#ifdef PR_SET_PTRACER //Required only for newer kernels that limit ptrace
		prctl(PR_SET_PTRACER, child_pid, 0, 0, 0);
		#endif
		close(fdPipe[0]);
		if(write(fdPipe[1], &message, 1) != 1)
		{	printf("Error communicating with debugger.\n");
			exit(code);
		}
		close(fdPipe[1]);
		waitpid(child_pid,NULL,0);
    }
    mpiUtil->exit(code);
}

// Stack trace for failed assertions
int assertStackTraceExit(const char* expr, const char* function, const char* file, long line)
{	fprintf(stderr, "%s:%ld: %s:\n\tAssertion '%s' failed", file, line, function, expr);
	gdbStackTraceExit(1);
	return 0;
}


namespace Citations
{
	//Single function whose static variable contains the list of citations
	//(Done this way, so that add may be called during static initialization safely)
	//If addCitation is non-null, enter the pair into the list
	//If getCitationList is non-null, retrieve the list
	void manage(std::pair<string,string>* addCitation=0, std::list<std::pair<string,string>>* getCitationList=0)
	{	static std::list<std::pair<string,string>> citationList; //pair.first = paper, pair.second = reason
		if(addCitation)
		{	auto iter=citationList.begin();
			bool foundPrev = false, duplicate = false;
			for(; iter!=citationList.end(); iter++)
			{	if(foundPrev && iter->second!=addCitation->second)
					break; //End of chain of citations with current paper
				if(iter->second==addCitation->second)
				{	foundPrev = true;
					if(iter->first==addCitation->first)
					{	duplicate = true;
						break;
					}
				}
			}
			if(!duplicate) citationList.insert(iter, *addCitation);
		}
		if(getCitationList) *getCitationList = citationList;
	}
	
	void add(string reason, string paper)
	{	std::pair<string,string> citation(paper, reason);
		manage(&citation, 0);
	}

	void print(FILE* fp)
	{	fprintf(fp, "\n---- Citations for features of the code used in this run ----\n\n");
		//Get the citation map:
		std::list<std::pair<string,string>> citationList;
		manage(0, &citationList);
		//Print entries:
		for(auto iter=citationList.begin(); iter!=citationList.end();)
		{	string paper = iter->first;
			for(;;iter++)
			{	if(iter==citationList.end() || iter->first != paper)
				{	//End the current chain of reasons for this paper
					istringstream iss(paper);
					while(!iss.eof())
					{	string line; getline(iss, line);
						fprintf(fp, "      %s\n", line.c_str());
					}
					fprintf(fp, "\n");
					break;
				}
				fprintf(fp, "   %s:\n", iter->second.c_str());
			}
		}
		fprintf(fp,
			"This list may not be complete. Please suggest additional citations and report\n"
			"any other bugs by creating a ticket at https://sourceforge.net/p/jdftx/tickets\n\n");
	}
}
