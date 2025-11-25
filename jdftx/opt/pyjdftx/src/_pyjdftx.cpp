#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#undef assert
#include <core/Util.h>
#include <electronic/run.h>
#include <electronic/JDFTxWrapper.h>


namespace py = pybind11;

extern int selectedDevice; //selected GPU number from core/GpuUtil.cpp (-1 if no GPU)

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

//View numpy arrays as equivalent NDarray object (to isolate libjdftx from python dependency)
NDarray viewFromPy(py::array_t<double> b)
{	py::buffer_info info = b.request();
	if(not (info.format == py::format_descriptor<double>::value))
		throw std::runtime_error("Incompatible buffer: need float64 (double) type");
	NDarray arr;
	arr.data = static_cast<double*>(info.ptr);
	for(size_t shape_i: info.shape) arr.shape.push_back(shape_i);
	for(size_t stride: info.strides) arr.strides.push_back(stride / sizeof(double));
	arr.onGpu = false;
	return arr;
}

//Reverse of above
py::array_t<double> viewToPy(const NDarray& arr)
{	std::string format = py::format_descriptor<double>::format();
	size_t n_dims = arr.shape.size();
	assert(not arr.onGpu);
	std::vector<py::ssize_t> shape;
	for(size_t shape_i: arr.shape) shape.push_back(shape_i);
	if(arr.strides.size())
	{	std::vector<py::ssize_t> byte_strides;
		for(size_t stride: arr.strides) byte_strides.push_back(sizeof(double) * stride);
		return py::array_t<double>(shape, byte_strides, arr.data, py::cast(arr)); //strided
	}
	else return py::array_t<double>(shape, arr.data, py::cast(arr)); //contiguous
}

#ifdef GPU_ENABLED
//Directly expose NDarray with __cuda_array_interface__ in binding for GPU code
typedef setCavity_t pySetCavity_t;
setCavity_t SetCavityWrapper(pySetCavity_t func) {return func;}
#else
//Wrap cavity callback from py::array to NDarray (to use np.ndarray in python CPU code)
typedef std::function<void(py::array_t<double>, std::vector<py::array_t<double>>)> pySetCavity_t;
struct SetCavityWrapper
{	pySetCavity_t func;
	SetCavityWrapper(pySetCavity_t func): func(func) {}
	void operator()(NDarray scalar, std::vector<NDarray> vec)
	{	std::vector<py::array_t<double>> pyVec;
		for(NDarray vec_i: vec) pyVec.push_back(viewToPy(vec_i));
		return func(viewToPy(scalar), pyVec);
	}
};
#endif

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
PYBIND11_MODULE(_pyjdftx_gpu, m) 
#else
PYBIND11_MODULE(_pyjdftx, m) 
#endif
{	if(import_mpi4py() < 0) throw std::runtime_error("Could not load mpi4py API.");
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
	
	m.def(
		"run", &run, py::arg("inputFilename"), py::arg("dryRun")=false, py::arg("printDefaults")=true,
		"Run a JDFTx calculation from an inputfile (just as the main executable would)."
	);
	
	m.def("selectedGPU", [](){return selectedDevice;}, "Selected GPU number, or -1 if CPU.");
	
	py::class_<NDarray>(m, "NDarray")
	#ifdef GPU_ENABLED
		.def_property_readonly("__cuda_array_interface__",
			[](NDarray& arr)
			{	assert(arr.onGpu);
				assert(not arr.strides.size()); //require contiguous for now
				py::dict result;		
				result["shape"] = py::tuple(py::cast(arr.shape));
				result["typestr"] = py::str(py::format_descriptor<double>::format());
				result["data"] = py::make_tuple(py::int_{reinterpret_cast<uintptr_t>(arr.data)}, false);
				result["version"] = 2;
				return result;
			}
		)
	#endif
	;
	
	py::class_<JDFTxWrapper>(m, "JDFTxWrapper")
		.def(
			py::init(
				[](std::vector<std::pair<string, string>> inputs, bool variableCell,
					std::pair<pySetCavity_t, pySetCavity_t> setCavity)
				{
					setCavity_t setCavityFuncs[2] = {setCavity_t(), setCavity_t()};
					if(setCavity.first)
					{	setCavityFuncs[0] = SetCavityWrapper(setCavity.first);
						setCavityFuncs[1] = SetCavityWrapper(setCavity.second);
					}
					return new JDFTxWrapper(inputs, variableCell, setCavityFuncs);
				}
			),
			py::arg("inputs"), py::arg("variableCell")=false,
			py::arg("setCavity")=std::make_pair(pySetCavity_t(), pySetCavity_t()),
			"__init__(self, inputs: list[tuple[str, str]], variableCell: bool)\n"
			"Setup calculation from jdftx command/value pairs in inputs."
		)
		.def("minimize", &JDFTxWrapper::minimize,
			 "Relax ions and/or lattice with parameters specified in commands."
		)
		.def("move", 
			 [](JDFTxWrapper& jw, py::array_t<double> delta_positions, py::array_t<double> delta_R)
			 {	jw.move(viewFromPy(delta_positions), viewFromPy(delta_R));
			 },
			 "move(self, delta_positions: np.ndarray, delta_R: np.ndarray)\n"
			 "Move ions and/or lattice by specified amounts.\n"
			 "Note that `delta_positions` is in fractional coordinates,\n"
			 "and `delta_R` is change of lattice vectors (in columns)."
		)
		.def("dumpEnd", &JDFTxWrapper::dumpEnd, "Dump final outputs at end of calculation")
		.def("getEnergy", &JDFTxWrapper::getEnergy, "Get current energy in Eh")
		.def("getForces",
			[](const JDFTxWrapper& jw) {return viewToPy(jw.getForces());}, 
			"Get current forces (Natoms x 3 array) in Eh/a0 (Cartesian)"
		)
		.def("getStress",
			[](const JDFTxWrapper& jw) {return viewToPy(jw.getStress());},
			"Get current stress (3 x 3 array) in Eh/a0^3 (Cartesian)"
		)
		.def_static("getCommands", &JDFTxWrapper::getCommands, "Get list of supported commands")
		.def_static("getCommandDoc", &JDFTxWrapper::getCommandDoc, "Get documentation of specified command");
}
