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
#include <getopt.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <commands/parser.h>

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#endif

#include <config.h> //This file is generated during build based on Git hash etc.

// Get the size of a file
off_t fileSize(const char *filename)
{	struct stat st;
	if(stat(filename, &st) == 0) return st.st_size;
    return -1;
}

//Program banner
void printVersionBanner()
{	logPrintf("\n*************** " PACKAGE_NAME " " VERSION_MAJOR_MINOR_PATCH
		" %s " PACKAGE_DESCRIPTION " ****************\n\n",
		(strlen(VERSION_HASH) ? "(git hash " VERSION_HASH ")" : ""));
}

//Print usage information
void printUsage(const char *name, const char* description)
{	printVersionBanner();
	logPrintf("Usage: %s [options]\n",name);
	logPrintf("\n\t%s\n\n", description);
	logPrintf("options:\n\n");
	logPrintf("\t-h --help               help (this output)\n");
	logPrintf("\t-v --version            version\n");
	logPrintf("\t-i --input <filename>   specify command input file, default = stdin\n");
	logPrintf("\t-o --output <filename>  specify output log file, default = stdout\n");
	logPrintf("\t-d --no-append          overwrite output file instead of appending\n");
	logPrintf("\t-t --template           print an input file template\n");
	logPrintf("\t-m --mpi-debug-log      write output from secondary MPI processes to jdftx.<proc>.mpiDebugLog (instead of /dev/null)\n");
	logPrintf("\t-n --dry-run            quit after initialization (to verify commands and other input files)\n");
	logPrintf("\t-c --cores              number of cores to use (ignored when launched using SLURM)\n");
	logPrintf("\t-s --skip-defaults      skip printing status of default commands issued automatically.\n");
	logPrintf("\n");
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
	stackTraceExit(1);
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
bool manualThreadCount = false;
static double startTime_us; //Time at which system was initialized in microseconds
const char* argv0 = 0;

void initSystem(int argc, char** argv)
{
	argv0 = argv[0]; //remember how the executable was issued (for stack traces)
	
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
	//---- hostname information
	std::vector<string> hostname(mpiUtil->nProcesses()); //list of hostnames by MPI process ID
	std::map< string, std::vector<int> > hostProcesses, hostGpuProcesses; //list of processes (and GPU-enabled processes) per hostname
	{	char hostnameTmp[256];
		gethostname(hostnameTmp, 256);
		hostname[mpiUtil->iProcess()] = hostnameTmp;
	}
	for(int jProcess=0; jProcess<mpiUtil->nProcesses(); jProcess++)
	{	mpiUtil->bcast(hostname[jProcess], jProcess);
		hostProcesses[hostname[jProcess]].push_back(jProcess);
		//Update GPU-enabled version of list:
		bool gpuEnabled = isGpuEnabled();
		mpiUtil->bcast(gpuEnabled, jProcess);
		if(gpuEnabled)
			hostGpuProcesses[hostname[jProcess]].push_back(jProcess);
	}
	logPrintf("Running on hosts (process indices):");
	for(const auto& iter: hostProcesses)
	{	logPrintf("  %s (", iter.first.c_str());
		for(int jProcess: iter.second) logPrintf(" %d", jProcess);
		logPrintf(" )");
	}
	logPrintf("\n");
	//---- commandline
	logPrintf("Executable %s with ", argv[0]);
	if(argc>1)
	{	logPrintf("command-line:");
		for(int i=1; i<argc; i++) logPrintf(" %s", argv[i]);
		logPrintf("\n");
	}
	else logPrintf("empty command-line (run with -h or --help for command-line options).\n");
	
	registerHandlers();
	
	double nGPUs = 0.;
	#ifdef GPU_ENABLED
	const std::vector<int>& gpuSiblings = hostGpuProcesses[hostname[mpiUtil->iProcess()]];
	if(!gpuInit(globalLog, &gpuSiblings, &nGPUs)) die_alone("gpuInit() failed\n\n")
	#endif
	
	//Divide up available cores between all MPI processes on a given node:
	if(!manualThreadCount) //skip if number of cores per process has been set with -c
	{	const std::vector<int>& siblings = hostProcesses[hostname[mpiUtil->iProcess()]];
		int nSiblings = siblings.size();
		int iSibling = std::find(siblings.begin(), siblings.end(), mpiUtil->iProcess()) - siblings.begin();
		nProcsAvailable = std::max(1, (nProcsAvailable * (iSibling+1))/nSiblings - (nProcsAvailable*iSibling)/nSiblings);
	}
	
	//Limit thread count if running within SLURM:
	const char* slurmCpusPerTask = getenv("SLURM_CPUS_PER_TASK");
	if(slurmCpusPerTask)
	{	int nThreadsMax;
		if(sscanf(slurmCpusPerTask, "%d", &nThreadsMax)==1)
			nProcsAvailable = nThreadsMax; //Slurm spec found, update available processor count (Thread.h)
		else
			logPrintf("Could not determine thread count from SLURM_CPUS_PER_TASK=\"%s\".\n", slurmCpusPerTask);
	}

	//Print number of threads per process:
	logPrintf("Maximum cpu threads by process:");
	for(int jProcess=0; jProcess<mpiUtil->nProcesses(); jProcess++)
	{	int nThreads = nProcsAvailable;
		mpiUtil->bcast(nThreads, jProcess);
		logPrintf(" %d", nThreads);
	}
	logPrintf("\n");
	resumeOperatorThreading(); //if necessary, this informs MKL of the thread count
	
	//Print total resources used by run:
	{	int nProcsTot = nProcsAvailable; mpiUtil->allReduce(nProcsTot, MPIUtil::ReduceSum);
		double nGPUsTot = nGPUs; mpiUtil->allReduce(nGPUsTot, MPIUtil::ReduceSum);
		logPrintf("Run totals: %d processes, %d threads, %lg GPUs\n", mpiUtil->nProcesses(), nProcsTot, nGPUsTot);
	}
	
	//Add citations to the code and general framework:
	Citations::add("Software package", "R. Sundararaman, D. Gunceler, K. Letchworth-Weaver, K.A. Schwarz and T.A. Arias, JDFTx, available from http://jdftx.sourceforge.net (2012)");
	Citations::add("Algebraic framework", "S. Ismail-Beigi and T.A. Arias, Computer Physics Communications 128, 1 (2000)");
}

void initSystemCmdline(int argc, char** argv, const char* description, string& inputFilename, bool& dryRun, bool& printDefaults)
{
	mpiUtil = new MPIUtil(argc, argv);
	
	//Parse command line:
	string logFilename; bool appendOutput=true;
	dryRun=false; printDefaults=true;
	option long_options[] =
		{	{"help", no_argument, 0, 'h'},
			{"version", no_argument, 0, 'v'},
			{"input",  required_argument, 0, 'i'},
			{"output", required_argument, 0, 'o'},
			{"no-append", no_argument, 0, 'd'},
			{"template", no_argument, 0, 't'},
			{"mpi-debug-log", no_argument, 0, 'm'},
			{"dry-run", no_argument, 0, 'n'},
			{"cores", required_argument, 0, 'c'},
			{"skip-defaults", no_argument, 0, 's'},
			{"write-manual", required_argument, 0, 'w'},
			{0, 0, 0, 0}
		};
	while (1)
	{	int c = getopt_long(argc, argv, "hvi:o:dtmnc:sw:", long_options, 0);
		if (c == -1) break; //end of options
		#define RUN_HEAD(code) if(mpiUtil->isHead()) { code } delete mpiUtil;
		switch (c)
		{	case 'v': RUN_HEAD( printVersionBanner(); ) exit(0);
			case 'h': RUN_HEAD( printUsage(argv[0], description); ) exit(0);
			case 'i': inputFilename.assign(optarg); break;
			case 'o': logFilename.assign(optarg); break;
			case 'd': appendOutput=false; break;
			case 't': RUN_HEAD( Everything e; printDefaultTemplate(e); ) exit(0);
			case 'm': mpiDebugLog=true; break;
			case 'n': dryRun=true; break;
			case 'c':
			{	int nCores = 0;
				if(sscanf(optarg, "%d", &nCores)==1 && nCores>0)
				{	nProcsAvailable=nCores;
					manualThreadCount =true;
				}
				break;
			}
			case 's': printDefaults=false; break;
			case 'w': RUN_HEAD( Everything e; writeCommandManual(e, optarg); ) exit(0);
			default: RUN_HEAD( printUsage(argv[0], description); ) exit(1);
		}
		#undef RUN_HEAD
	}
	
	// Open the logfile (if any):
	if(logFilename.length())
	{	globalLog = fopen(logFilename.c_str(), appendOutput ? "a" : "w");
		if(!globalLog)
		{	globalLog = stdout;
			logPrintf("WARNING: Could not open log file '%s' for writing, using standard output.\n", logFilename.c_str());
		}
	}

	//Print banners, setup threads, GPUs and signal handlers
	initSystem(argc, argv);
}

#ifdef ENABLE_PROFILING
void stopWatchManager(const StopWatch* addWatch=0, const string* watchName=0)
{	static std::multimap<string, const StopWatch*> watches; //static array of all watches
	if(addWatch) watches.insert(std::make_pair(*watchName,addWatch));
	else //print timings:
	{	logPrintf("\n");
		for(const auto& wPair: watches) wPair.second->print();
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
	logPrintf("\n");
	ManagedMemory::reportUsage();
	#endif
	
	if(!mpiUtil->isHead())
	{	if(mpiDebugLog) fclose(globalLog);
		globalLog = 0;
	}
	fclose(nullLog);
	if(globalLog && globalLog != stdout)
		fclose(globalLog);
	delete mpiUtil;
}


//------------ Timing helpers ----------------

inline double clock_us_epoch() //time in microseconds (since standard POSIX specified epoch)
{	timeval tv;
	gettimeofday(&tv,NULL);
	return ((tv.tv_sec & 0x1fffff) * 1e6) + tv.tv_usec;
}
double clock_us()
{	static double tStart = clock_us_epoch();
	return clock_us_epoch() - tStart;
}
double clock_sec()
{	return 1e-6*clock_us();
}


#ifdef ENABLE_PROFILING
StopWatch::StopWatch(string name) : Ttot(0), TsqTot(0), nT(0), name(name) { stopWatchManager(this, &name); }
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
void printStack(bool detailedStackScript)
{	const int maxStackLength = 1024;
	void* tracePtrs[maxStackLength];
	int count = backtrace(tracePtrs, maxStackLength);
	char** funcNames = backtrace_symbols(tracePtrs, count);
	logPrintf("\nStack trace:\n");
	for(int i=0; i<count; i++)
		logPrintf("\t%2d: %s\n", i, funcNames[i]);
	if(detailedStackScript)
	{	logPrintf("Writing 'jdftx-stacktrace' (for use with script printStackTrace): "); logFlush();
		FILE* fp = fopen("jdftx-stacktrace", "w");
		if(fp)
		{	for(int i=0; i<count; i++)
				fprintf(fp, "%s\n", funcNames[i]);
			fclose(fp);
			logPrintf("done.\n");
		}
		else logPrintf("could not open file for writing.\n");
	}
	free(funcNames);
}

// Exit on error with a more in-depth stack trace
void stackTraceExit(int code)
{	printStack(true);
	mpiUtil->exit(code);
}

// Stack trace for failed assertions
int assertStackTraceExit(const char* expr, const char* function, const char* file, long line)
{	fprintf(stderr, "%s:%ld: %s:\n\tAssertion '%s' failed", file, line, function, expr);
	stackTraceExit(1);
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
