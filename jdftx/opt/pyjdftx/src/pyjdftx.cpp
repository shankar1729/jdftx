#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#undef assert
#include <core/Util.h>
#include <electronic/JDFTxWrapper.h>


namespace py = pybind11;

// Wrapper to convert between MPI_Comm and mpi4pyComm
struct mpi4pyComm
{	mpi4pyComm() = default;
	mpi4pyComm(MPI_Comm comm) : comm(comm) {}
	operator MPI_Comm() {return comm;}
	MPI_Comm comm;
};

// Register wrapper with pybind11
namespace pybind11 { namespace detail {
	template <> struct type_caster<mpi4pyComm>
	{	PYBIND11_TYPE_CASTER(mpi4pyComm, _("mpi4pyComm"));
		bool load(handle src, bool) 
		{	PyObject *py_src = src.ptr();
			if(PyObject_TypeCheck(py_src, &PyMPIComm_Type)) 
				value.comm = *PyMPIComm_Get(py_src);
			else return false;
			return !PyErr_Occurred();
		}
		static handle cast(mpi4pyComm src, return_value_policy, handle)
		{	return PyMPIComm_New(src.comm);
		}
  };
}}

void initialize(mpi4pyComm comm, mpi4pyComm commAll, string logFilename, bool appendLog)
{	mpiWorldFull = new MPIUtil(commAll);
	mpiWorld = new MPIUtil(comm); //done later, so that iProc from here fixes random seed
	globalLog = fopen(logFilename.c_str(), appendLog ? "a" : "w");
	if(!globalLog)
	{	globalLog = stdout;
		logPrintf("WARNING: Could not open log file '%s' for writing, using standard output.\n", logFilename.c_str());
	}
	int argc = 1;
	char* argv[] = {(char*)JDFTX_LIBRARY};  //reported in log and stack traces
	initSystem(argc, argv);
}


#ifdef GPU_ENABLED
PYBIND11_MODULE(pyjdftx_gpu, m) 
#else
PYBIND11_MODULE(pyjdftx, m) 
#endif
{	if (import_mpi4py() < 0) {
		throw std::runtime_error("Could not load mpi4py API.");
	}
	m.doc() = "Python wrapper to JDFTx";
	m.def(
		"initialize", &initialize,
		"initialize(comm: MPI.Comm, commAll: MPI.Comm, logFilename:str, appendLog: bool)\n"
		"Initialize hardware resources including MPI and GPUs (if any).\n"
		"Here, `comm` is the communicator to use during the run,\n"
		"while `commAll` is the overall communicator over which initialize\n"
		"is being called simultaneously in order to divide cores/GPUs correctly.\n\n"
		"If commAll is set to comm during a split run, it is the responsibility\n"
		"of the calling script to set the correct number of CPU threads to use\n"
		"in SLURM_CPUS_PER_TASK, and to modify CUDA_VISIBLE_DEVICES to select\n"
		"the GPU that each process should use (or all that the current comm\n"
		"should have access to, in order to avoid overcommitting resources."
	);
	
	m.def(
		"finalize", &finalizeSystem,
		"finalize(success: bool)\n"
		"Report resource usage and clean up MPI, logs etc."
	);
	
	py::class_<JDFTxWrapper>(m, "JDFTxWrapper")
        .def(
			py::init<std::vector<std::pair<string, string>>, bool>(),
			"__init__(self, inputs: list[tuple[str, str]], variableCell: bool)\n"
			"Setup calculation from jdftx command/value pairs in inputs."
		)
		.def("minimize", &JDFTxWrapper::minimize,
			 "Relax ions and/or lattice with parameters specified in commands."
		 );
}
