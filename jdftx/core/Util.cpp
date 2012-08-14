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
#include <sys/prctl.h>

#ifdef GPU_ENABLED
#include <core/GpuUtil.h>
#endif

#include <config.h> //This file is generated during build based on SVN version etc.

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
{	resetHandlers();
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
		{	case 'q': case 'Q': printf("Quitting now ...\n"); exit(0);
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

void initSystem(int argc, char** argv)
{
	//Print a welcome banner with useful information
	printVersionBanner();
	time_t timenow = time(0);
	logPrintf("Start date and time: %s", ctime(&timenow)); //note ctime output has a "\n" at the end
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
	const char* slurmCpusPerNode = getenv("SLURM_JOB_CPUS_PER_NODE");
	if(slurmCpusPerNode)
	{	int nThreadsMax;
		if(sscanf(slurmCpusPerNode, "%d", &nThreadsMax)==1)
			nProcsAvailable = nThreadsMax; //Slurm spec found, update available processor count (Thread.h)
		else
			logPrintf("Could not determine thread count from SLURM_JOB_CPUS_PER_NODE=\"%s\".\n", slurmCpusPerNode);
	}
	logPrintf("Will run with a maximum of %d cpu threads.\n", nProcsAvailable);
}



#ifdef ENABLE_PROFILING
StopWatch::StopWatch(string name) : Ttot(0), TsqTot(0), nT(0), name(name) {}
void StopWatch::start() { tPrev = clock_us(); }
void StopWatch::stop() { double T = clock_us()-tPrev; Ttot+=T; TsqTot+=T*T; nT++; }
StopWatch::~StopWatch()
{	if(nT)
	{	double meanT = Ttot/nT;
		double sigmaT = sqrt(TsqTot/nT - meanT*meanT);
		printf("'%s' took %lf +/- %lf s (called %d times for a total of %lf s)\n",
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
	int fdPipe[2]; if(pipe(fdPipe)) { printf("Error creating pipe.\n"); exit(code); }
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
	{	prctl(PR_SET_PTRACER, child_pid, 0, 0, 0);
		close(fdPipe[0]);
		if(write(fdPipe[1], &message, 1) != 1)
		{	printf("Error communicating with debugger.\n");
			exit(code);
		}
		close(fdPipe[1]);
		waitpid(child_pid,NULL,0);
    }
    exit(code);
}

// Stack trace for failed assertions
int assertStackTraceExit(const char* expr, const char* function, const char* file, long line)
{	fprintf(stderr, "%s:%ld: %s:\n\tAssertion '%s' failed", file, line, function, expr);
	gdbStackTraceExit(1);
	return 0;
}


FILE* globalLog = stdout; // this might be replaced by a file pointer in main

