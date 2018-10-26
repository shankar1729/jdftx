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
#include <core/ManagedMemory.h>
#include <core/GpuUtil.h>
#include <cmath>
#include <csignal>
#include <list>
#include <algorithm>
#include <getopt.h>
#include <commands/parser.h>

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#endif

#include <config.h> //This file is generated during build based on Git hash etc.

InitParams::InitParams(const char* description, class Everything* e)
: description(description), e(e), packageName(0), versionString(0), versionHash(0)
{
}

//Program banner
void printVersionBanner(const InitParams* ip)
{	string deco(15, '*');
	string prefix;
	logPrintf("\n");
	if(ip && ip->packageName)
	{	//Print other package's information (JDFTx is only a helper)
		logPrintf("%s %s %s %s %s\n", deco.c_str(), ip->packageName, ip->versionString,
			(strlen(ip->versionHash) ? ("(git hash " + string(ip->versionHash) + ")").c_str() : ""), deco.c_str());
		deco.assign(15, '+');
		prefix = "Linked to ";
	}
	logPrintf("%s %s" PACKAGE_NAME " " VERSION_MAJOR_MINOR_PATCH " %s %s\n",
		deco.c_str(), prefix.c_str(), (strlen(VERSION_HASH) ? "(git hash " VERSION_HASH ")" : ""), deco.c_str());
	logPrintf("\n"); logFlush();
}

//Print usage information
void printUsage(const char *name, const InitParams& ip)
{	printVersionBanner(&ip);
	logPrintf("Usage: %s [options]\n",name);
	logPrintf("\n\t%s\n\n", ip.description);
	logPrintf("options:\n\n");
	logPrintf("\t-h --help               help (this output)\n");
	logPrintf("\t-v --version            version\n");
	logPrintf("\t-i --input <filename>   specify command input file, default = stdin\n");
	logPrintf("\t-o --output <filename>  specify output log file, default = stdout\n");
	logPrintf("\t-d --no-append          overwrite output file instead of appending\n");
	logPrintf("\t-t --template           print an input file template\n");
	logPrintf("\t-m --mpi-debug-log      write output from secondary MPI processes to jdftx.<proc>.mpiDebugLog (instead of /dev/null)\n");
	logPrintf("\t-n --dry-run            quit after initialization (to verify commands and other input files)\n");
	logPrintf("\t-c --cores              number of cores per process (ignored when launched using SLURM)\n");
	logPrintf("\t-G --nGroups            number of MPI process groups (default or 0 => each process in own group of size 1)\n");
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
{	if(feof(stdin)) mpiWorld->exit(0);
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
		{	case 'q': case 'Q': printf("Quitting now ...\n"); mpiWorld->exit(0);
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

string inputBasename("stdin"); //!< Basename of input file or "stdin" that can be used as a default run-name
FILE* globalLog = stdout; // this might be replaced by a file pointer in main before calling initSystem()
FILE* globalLogOrig;
FILE* nullLog = 0;

void logSuspend()
{	if(nullLog) globalLog = nullLog;
}

void logResume()
{	globalLog = globalLogOrig;
}

int nProcessGroups = 0;
MPIUtil* mpiWorld = 0;
MPIUtil* mpiGroup = 0;
MPIUtil* mpiGroupHead = 0;
bool mpiDebugLog = false;
bool manualThreadCount = false;
size_t mempoolSize = 0;
static double startTime_us; //Time at which system was initialized in microseconds
const char* argv0 = 0;
uint32_t crc32(const string& s); //CRC32 checksum for a string (implemented below)

//Print process index distribution given communicators:
void printProcessDistribution(string header, string label, const MPIUtil* mpiUtil, const MPIUtil* mpiUtilHead)
{	if(mpiUtil->isHead())
	{	ostringstream oss;
		oss << label << " (";
		int jStart=mpiWorld->iProcess(), jStop=jStart; //range yet to be printed
		for(int jProc=1; jProc<mpiUtil->nProcesses(); jProc++)
		{	int jCur; //rank of jProc in world
			mpiUtil->recv(jCur, jProc, 0);
			if(jCur == jStop+1)
				jStop = jCur; //contiguous with previous range
			else
			{	oss << jStart; if(jStop>jStart) oss << '-' << jStop; oss << ','; //flush previous range
				jStart = jStop = jCur; //start new range
			}
		}
		oss << jStart; if(jStop>jStart) oss << '-' << jStop; oss << ')'; //flush pending range
		//send to world head to print:
		if(mpiUtilHead->isHead()) //note that mpiUtilHead->head == mpiWorld->head
		{	logPrintf("%s (process indices):  %s", header.c_str(), oss.str().c_str());
			for(int jHead=1; jHead<mpiUtilHead->nProcesses(); jHead++)
			{	string buf;
				mpiUtilHead->recv(buf, jHead, 0);
				logPrintf("  %s", buf.c_str());
			}
			logPrintf("\n");
		}
		else mpiUtilHead->send(oss.str(), 0, 0);
	}
	else mpiUtil->send(mpiWorld->iProcess(), 0, 0);
}

void initSystem(int argc, char** argv, const InitParams* ip)
{
	argv0 = argv[0]; //remember how the executable was issued (for stack traces)
	
	if(!mpiWorld) mpiWorld = new MPIUtil(argc, argv);
	nullLog = fopen("/dev/null", "w");
	if(!mpiWorld->isHead())
	{	if(mpiDebugLog)
		{	char fname[256]; sprintf(fname, "jdftx.%d.mpiDebugLog", mpiWorld->iProcess());
			globalLog = fopen(fname, "w");
		}
		else globalLog = nullLog;
	}
	globalLogOrig = globalLog;
	
	//Star time and commandline:
	printVersionBanner(ip);
	time_t startTime = time(0);
	startTime_us = clock_us();
	logPrintf("Start date and time: %s", ctime(&startTime)); //note ctime output has a "\n" at the end
	//---- commandline
	logPrintf("Executable %s with ", argv[0]);
	if(argc>1)
	{	logPrintf("command-line:");
		for(int i=1; i<argc; i++) logPrintf(" %s", argv[i]);
		logPrintf("\n");
	}
	else logPrintf("empty command-line (run with -h or --help for command-line options).\n");
	registerHandlers();

	//Determine and print distribution of processes on hosts:
	//--- get current hostname and checksum:
	string hostname;
	{	char hostnameTmp[256];
		gethostname(hostnameTmp, 256);
		hostname = string(hostnameTmp);
	}
	int hostsum = abs(int(crc32(hostname))); //ensure positive for MPI_split below
	//---- create MPI group for each host:
	MPIUtil mpiHost(0,0, MPIUtil::ProcDivision(mpiWorld, 0, hostsum));
	MPIUtil mpiHostGpu(0,0, MPIUtil::ProcDivision(mpiWorld, 0, isGpuEnabled() ? hostsum : 0)); //for grouping processes with GPUs
	MPIUtil mpiHostHead(0,0, MPIUtil::ProcDivision(mpiWorld, 0, mpiHost.iProcess())); //communicator between similar rank within each host
	printProcessDistribution("Running on hosts", hostname, &mpiHost, &mpiHostHead);
	
	//Initialize process groups:
	if(nProcessGroups <= 0) nProcessGroups = mpiWorld->nProcesses(); //default: one group per process
	mpiGroup = new MPIUtil(0,0, MPIUtil::ProcDivision(mpiWorld, nProcessGroups));
	mpiGroupHead = new MPIUtil(0,0, MPIUtil::ProcDivision(mpiWorld, 0, mpiGroup->iProcess())); //communicator between similar rank within each group
	{	ostringstream oss; oss << mpiGroup->procDivision.iGroup;
		printProcessDistribution("Divided in process groups", oss.str(), mpiGroup, mpiGroupHead);
	}
	
	double nGPUs = 0.;
	#ifdef GPU_ENABLED
	if(!gpuInit(globalLog, &mpiHostGpu, &nGPUs)) die_alone("gpuInit() failed\n\n")
	#endif
	
	//Override available cores per node if specified:
	const char* envCpusPerNode = getenv("JDFTX_CPUS_PER_NODE");
	if(envCpusPerNode)
	{	int nCpusPerNode;
		if(sscanf(envCpusPerNode, "%d", &nCpusPerNode)==1)
			nProcsAvailable = nCpusPerNode; //Override found, update available processor count (Thread.h)
		else
			logPrintf("Could not determine total core count from JDFTX_CPUS_PER_NODE=\"%s\".\n", envCpusPerNode);
	}
	
	//Divide up available cores between all MPI processes on a given node:
	if(!manualThreadCount) //skip if number of cores per process has been set with -c
	{	int nSiblings = mpiHost.nProcesses();
		int iSibling = mpiHost.iProcess();
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
	resumeOperatorThreading(); //if necessary, this informs MKL of the thread count
	
	//Print total resources used by run:
	{	int nProcsTot = nProcsAvailable; mpiWorld->allReduce(nProcsTot, MPIUtil::ReduceSum);
		double nGPUsTot = nGPUs; mpiWorld->allReduce(nGPUsTot, MPIUtil::ReduceSum);
		logPrintf("Resource initialization completed at t[s]: %9.2lf\n", clock_sec());
		logPrintf("Run totals: %d processes, %d threads, %lg GPUs\n", mpiWorld->nProcesses(), nProcsTot, nGPUsTot);
	}
	
	//Memory pool size:
	const char* mempoolSizeStr = getenv("JDFTX_MEMPOOL_SIZE");
	if(mempoolSizeStr)
	{	int mempoolSizeMB;
		if(sscanf(mempoolSizeStr, "%d", &mempoolSizeMB)==1 && mempoolSizeMB>=0)
		{	mempoolSize = ((size_t)mempoolSizeMB) << 20; //convert to bytes
			logPrintf("Memory pool size: %d MB (per process)\n", mempoolSizeMB);
		}
		else
			logPrintf("Could not determine memory pool size from JDFTX_MEMPOOL_SIZE=\"%s\".\n", mempoolSizeStr);
	}
	
	//Add citations to the code for all calculations:
	Citations::add("Software package",
		"R. Sundararaman, K. Letchworth-Weaver, K.A. Schwarz, D. Gunceler, Y. Ozhabes and T.A. Arias, "
		"'JDFTx: software for joint density-functional theory', SoftwareX 6, 278 (2017)");
}

void initSystemCmdline(int argc, char** argv, InitParams& ip)
{
	mpiWorld = new MPIUtil(argc, argv);
	
	//Parse command line:
	string logFilename; bool appendOutput=true;
	ip.dryRun=false; ip.printDefaults=true;
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
			{"nGroups", required_argument, 0, 'G'},
			{"skip-defaults", no_argument, 0, 's'},
			{"write-manual", required_argument, 0, 'w'},
			{0, 0, 0, 0}
		};
	while (1)
	{	int c = getopt_long(argc, argv, "hvi:o:dtmnc:G:sw:", long_options, 0);
		if (c == -1) break; //end of options
		#define RUN_HEAD(code) if(mpiWorld->isHead()) { code } delete mpiWorld;
		switch (c)
		{	case 'v': RUN_HEAD( printVersionBanner(&ip); ) exit(0);
			case 'h': RUN_HEAD( printUsage(argv[0], ip); ) exit(0);
			case 'i': ip.inputFilename.assign(optarg); break;
			case 'o': logFilename.assign(optarg); break;
			case 'd': appendOutput=false; break;
			case 't': RUN_HEAD( if(ip.e) printDefaultTemplate(*ip.e); ) exit(0);
			case 'm': mpiDebugLog=true; break;
			case 'n': ip.dryRun=true; break;
			case 'c':
			{	int nCores = 0;
				if(sscanf(optarg, "%d", &nCores)==1 && nCores>0)
				{	nProcsAvailable=nCores;
					manualThreadCount =true;
				}
				break;
			}
			case 'G':
			{	if(!(sscanf(optarg, "%d", &nProcessGroups)==1 && nProcessGroups>=0))
				{	RUN_HEAD(
						printf("\nOption -G (--nGroups) must be a non-negative integer.\n");
						printUsage(argv[0], ip);
					)
					exit(1);
				}
				break;
			}
			case 's': ip.printDefaults=false; break;
			case 'w': RUN_HEAD( if(ip.e) writeCommandManual(*ip.e, optarg); ) exit(0);
			default: RUN_HEAD( printUsage(argv[0], ip); ) exit(1);
		}
		#undef RUN_HEAD
	}
	
	//Open the logfile (if any):
	if(logFilename.length())
	{	globalLog = fopen(logFilename.c_str(), appendOutput ? "a" : "w");
		if(!globalLog)
		{	globalLog = stdout;
			logPrintf("WARNING: Could not open log file '%s' for writing, using standard output.\n", logFilename.c_str());
		}
	}
	
	//Set input base name if necessary:
	if(ip.inputFilename.length())
	{	inputBasename = ip.inputFilename;
		//Remove extension:
		size_t lastDot = inputBasename.find_last_of(".");
		if(lastDot != string::npos)
			inputBasename = inputBasename.substr(0, lastDot); //Remove extension
		//Remove leading path:
		size_t lastSlash = inputBasename.find_last_of("\\/");
		if(lastSlash != string::npos)
			inputBasename = inputBasename.substr(lastSlash+1);
	}
	
	//Print banners, setup threads, GPUs and signal handlers
	initSystem(argc, argv, &ip);
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
		if(mpiWorld->isHead() && globalLog != stdout)
			fprintf(stderr, "Failed.\n");
	}
	
	#ifdef ENABLE_PROFILING
	stopWatchManager();
	logPrintf("\n");
	ManagedMemoryBase::reportUsage();
	#endif
	
	if(!mpiWorld->isHead())
	{	if(mpiDebugLog) fclose(globalLog);
		globalLog = 0;
	}
	fclose(nullLog);
	if(globalLog && globalLog != stdout)
		fclose(globalLog);
	delete mpiGroupHead;
	delete mpiGroup;
	delete mpiWorld;
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

//--------------- Miscellaneous ------------------

// Get the size of a file
off_t fileSize(const char *filename)
{	struct stat st;
	if(stat(filename, &st) == 0) return st.st_size;
    return -1;
}

bool isLittleEndian()
{	static bool isLE = false, initializedLE = false;
	if(!initializedLE)
	{	//To be absolutely sure, check with an 8-byte type:
		uint64_t testData = 0x0001020304050607;
		uint8_t* testPtr = (uint8_t*)(&testData);
		bool littleCheck = true;
		bool bigCheck = true;
		for(int j=0; j<8; j++)
		{	littleCheck &= (testPtr[j]==7-j);
			bigCheck &= (testPtr[j]==j);
		}
		if(littleCheck) isLE = true; //little endian
		else if(bigCheck) isLE = false; //big endian
		else //neither:
			die("Binary I/O not yet supported on mixed-endian CPUs.\n")
	}
	return isLE;
}

//Convert data from operating endianness to little-endian
void convertToLE(void* ptr, size_t size, size_t nmemb)
{	static StopWatch watch("endianSwap");
	if(isLittleEndian() || size==1) return; //nothing to do on little-endian systems, or on byte-wise data
	watch.start();
	//Determine chunk size over which byte-swapping must occur:
	size_t chunkSize = 0, nChunks = 0;
	switch(size)
	{	case 2:
		case 4:
		case 8:
			chunkSize = size;
			nChunks = nmemb;
			break;
		default:
			if(size % 8 == 0)
			{	//assume composite data set of doubles eg. complex
				chunkSize = 8;
				nChunks = nmemb * (size/8);
			}
			else
				die("Unsupported size '%zu' for binary I/O on big-endian systems.\n", size)
	}
	//Apply byte-swapping:
	char* bytes = (char*)ptr;
	for(size_t iChunk=0; iChunk<nChunks; iChunk++)
	{	for(size_t iByte=0; iByte<chunkSize/2; iByte++)
		{	std::swap(bytes[iByte], bytes[chunkSize-1-iByte]);
		}
		bytes += chunkSize; //move to next chunk
	}
	watch.stop();
}

//Convert data from little-endian to operating endianness
void convertFromLE(void* ptr, size_t size, size_t nmemb)
{	convertToLE(ptr, size, nmemb); //swapping between little and big endian is its own inverse
}

//Read from a little-endian binary file, regardless of operating endianness
size_t freadLE(void *ptr, size_t size, size_t nmemb, FILE* fp)
{	size_t result = fread(ptr, size, nmemb, fp); //byte-by-byte read
	convertFromLE(ptr, size, nmemb); //modify byte-order in place, if necessary
	return result;
}

//Write to a little-endian binary file, regardless of operating endianness
size_t fwriteLE(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{	convertToLE((void*)ptr, size, nmemb); //modify byte-order in place, if necessary
	size_t result = fwrite((void*)ptr, size, nmemb, fp); //byte-by-byte write
	convertFromLE((void*)ptr, size, nmemb); //restore byte-order (since data should not be logically modified)
	return result;
}


// Exit on error with a more in-depth stack trace
void stackTraceExit(int code)
{	printStack(true);
	mpiWorld->exit(code);
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
			"This list may not be complete. Please suggest additional citations or\n"
			"report any other bugs at https://github.com/shankar1729/jdftx/issues\n\n");
		fflush(fp);
	}
}


//Calculate CRC32 checksum for a buffer
template<typename CharIterator> uint32_t crc32(const CharIterator begin, const CharIterator end)
{	//--- Code based on CRC32 snippet from c.snippets.org ---
	static uint32_t table[] = {
		0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3, 0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
		0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91, 0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
		0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5, 0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
		0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b, 0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
		0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
		0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d, 0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
		0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01, 0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
		0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
		0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
		0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
		0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
		0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683, 0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
		0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
		0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5, 0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
		0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
		0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
		0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
		0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
		0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45, 0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
		0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
		0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
		0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d };
	//--- calculate and return checksum:
	uint32_t result = 0xFFFFFFFF;
	for(CharIterator ptr=begin; ptr<end; ptr++)
		result = (table[(result^(*ptr)) & 0xff] ^ (result >> 8));
	return ~result;
}
uint32_t crc32(const string& s)
{	return crc32(s.begin(), s.end());
}
