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

void initialize(mpi4pyComm comm, string logFilename, bool appendLog=true)
{	mpiWorld = new MPIUtil(comm);
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
	m.doc() = "pybind11 example plugin"; // optional module docstring
	m.def("initialize", &initialize, "Initialize hardware resources including MPI and GPUs (if any)");
	
	py::class_<JDFTxWrapper>(m, "JDFTxWrapper")
        .def(py::init<std::vector<std::pair<string, string>>>());
}
