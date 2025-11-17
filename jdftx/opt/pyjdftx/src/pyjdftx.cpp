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
	template<> struct type_caster<mpi4pyComm>
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

//Convert python buffers to equivalent NDarray object (to isolate libjdftx from python dependency)
NDarray convert(py::buffer b)
{	py::buffer_info info = b.request();
	if(not (info.format == py::format_descriptor<double>::value))
		throw std::runtime_error("Incompatible buffer: need float64 (double) type");
	NDarray arr;
	arr.data = static_cast<double*>(info.ptr);
	for(size_t shape_i: info.shape) arr.shape.push_back(shape_i);
	for(size_t stride: info.strides) arr.strides.push_back(stride / sizeof(double));
	return arr;
}

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
		)
		.def("move", 
			 [](JDFTxWrapper& jw, py::buffer delta_positions, py::buffer delta_R)
			 {	jw.move(convert(delta_positions), convert(delta_R));
			 },
			 "move(self, delta_positions: np.ndarray, delta_R: np.ndarray)\n"
			 "Move ions and/or lattice by specified amounts.\n"
			 "Note that `delta_positions` is in fractional coordinates,\n"
			 "and `delta_R` is change of lattice vectors (in columns)."
		)
		.def("getEnergy", &JDFTxWrapper::getEnergy, "Get current energy in Eh")
		.def("getForces", &JDFTxWrapper::getForces, "Get current forces in Eh/a0 (Cartesian)")
		.def("getStress", &JDFTxWrapper::getStress, "Get current stress in Eh/a0^3 (Cartesian)");
	
	py::class_<NDarray>(m, "NDarray", py::buffer_protocol())
		.def_buffer(
			[](NDarray& arr) -> py::buffer_info 
			{
				std::string format = py::format_descriptor<double>::format();
				size_t n_dims = arr.shape.size();
				std::vector<py::ssize_t> shape, byte_strides;
				for(size_t shape_i: arr.shape) shape.push_back(shape_i);
				for(size_t stride: arr.strides) byte_strides.push_back(sizeof(double) * stride);
				return py::buffer_info((double*)arr.data, //cast to non-const pointer
					sizeof(double), format, n_dims, shape, byte_strides, true); //but marked read-only
			}
		);
}
