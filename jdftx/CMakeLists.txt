cmake_minimum_required(VERSION 2.8.4)

project(JDFTx)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake-Modules/")

#Package configuration:
set(CPACK_PACKAGE_NAME "JDFTx")
set(CPACK_PACKAGE_DESCRIPTION "The playground for joint density functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Toolkit for rapid development of density-functional methods")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "5")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

#Get revision information
find_package(Git)
if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE IS_GIT
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if("true" STREQUAL "${IS_GIT}")	
		execute_process(
			COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_HASH
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		message(STATUS "Git revision hash: ${GIT_HASH}")
	endif()
endif()

#Static linking:
option(StaticLinking "Link statically instead of dynamically" OFF)
if(CYGWIN) #Force static linking on Windows:
	if(NOT StaticLinking)
		message(STATUS "Forcing static linking for Cygwin/Windows")
		set(StaticLinking ON)
	endif()
endif()
if(StaticLinking)
	if(WIN32)
		set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	else()
		set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif()
	set(CMAKE_EXE_LINKER_FLAGS -static)
	set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
	set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
	set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
	set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
	set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
	set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif()

#GSL:
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

#Use indepedent CBLAS, LAPACK and FFTW (default) or MKL for all three
option(EnableMKL "Use Intel MKL to provide BLAS, LAPACK and FFTs")
option(ForceFFTW "Force usage of FFTW (even if MKL is enabled)")
option(ThreadedBLAS "Used built-in threading of the BLAS library if yes; thread in JDFTx if no (currently affects only MKL)" ON)
option(EnableScaLAPACK "Enable ScaLAPACK support (currently used only by the BerkeleyGW output option)")
set(CMAKE_THREAD_PREFER_PTHREAD)
find_package(Threads REQUIRED)
if(EnableMKL)
	find_package(MKL REQUIRED)
	set(CBLAS_LAPACK_FFT_LIBRARIES ${MKL_LIBRARIES}) #MKL provides CBLAS, FFTW3 and LAPACK
	add_definitions("-DMKL_PROVIDES_BLAS") #Special handling is required for BLAS calls from threaded sections
	include_directories(${MKL_INCLUDE_DIR})
		if(ForceFFTW)
		find_package(FFTW3 REQUIRED)
		set(CBLAS_LAPACK_FFT_LIBRARIES ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${MKL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) #Explicit FFTW3, rest from MKL
	else()
		add_definitions("-DMKL_PROVIDES_FFT") #Special handling is required for FFT initialization
		set(CBLAS_LAPACK_FFT_LIBRARIES ${MKL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) #MKL provides CBLAS, FFTW3 and LAPACK
	endif()
	if(ThreadedBLAS)
		add_definitions("-DTHREADED_BLAS")
	endif()
else()
	find_package(FFTW3 REQUIRED)
	find_package(LAPACK_ATLAS REQUIRED)
	find_package(CBLAS REQUIRED)
	set(CBLAS_LAPACK_FFT_LIBRARIES ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${CBLAS_LIBRARY} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()
include_directories(${FFTW3_INCLUDE_DIR})

option(EnableMPI "Use MPI parallelization (in addition to threads / gpu)" ON)
if(EnableMPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH};${MPI_INCLUDE_PATH}") #To find libraries installed to MPI specific paths
	add_definitions("-DMPI_ENABLED")
endif()

option(EnableLibXC "Use LibXC to provide additional exchange-correlation functionals")
if(EnableLibXC)
	find_package(LIBXC REQUIRED)
	include_directories(${LIBXC_INCLUDE_DIR})
	add_definitions("-DLIBXC_ENABLED")
	#Generate functional list available from LibXC:
	execute_process(COMMAND ${CMAKE_SOURCE_DIR}/opt/indexLibXC.sh ${LIBXC_INCLUDE_DIR}/xc_funcs.h OUTPUT_FILE ${CMAKE_BINARY_DIR}/xcMap.h)
endif()

option(EnableHDF5 "Enable HDF5 features (currently required only by Berkeley GW dump option)")
if(EnableHDF5)
	find_package(HDF5 REQUIRED)
	include_directories(${HDF5_INCLUDE_DIRS})
	add_definitions("-DHDF5_ENABLED")
endif()

#Process configuration information into config.h (with config.in.h as a template)
configure_file(${CMAKE_SOURCE_DIR}/config.in.h ${CMAKE_BINARY_DIR}/config.h)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR})

option(EnableProfiling "Optionally enable profiling to determine ideal functions to optimize (for developers)")
if(EnableProfiling)
	add_definitions("-DENABLE_PROFILING")
endif()

#------------------------- Compilation flags ------------------------

#Base flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -Wall -O3")
set(JDFTX_CPU_FLAGS "-std=c++0x") #Flags used for host compiler alone (i.e. not for gpu device code)

#Additional compiler-dependent flags:
include(CheckCXXCompilerFlag)
#--- Suppress gcc warning
check_cxx_compiler_flag(-Wno-unused-result HAS_NO_UNUSED_RESULT)
if(HAS_NO_UNUSED_RESULT)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
endif()
option(LinkTimeOptimization "Enable link-time optimizations (eg. flto on GNU, ipo on Intel).")
if(LinkTimeOptimization)
	#--- Link-time optimizations (GCC / CLang)
	if(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
		check_cxx_compiler_flag(-flto HAS_FLTO)
		if(HAS_FLTO)
			set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -flto")
		endif()
	endif()
	#--- Link-time optimizations (ICC)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		check_cxx_compiler_flag(-ipo HAS_IPO)
		if(HAS_IPO)
			set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -ipo")
		endif()
	endif()
endif()
#--- Template recursion depth (more needed in Clang)
check_cxx_compiler_flag(-ftemplate-depth-512 HAS_TEMPLATE_DEPTH)
if(HAS_TEMPLATE_DEPTH)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-512")
endif()

#Architecture dependent optimizations
option(CompileNative "Enable aggressive architecture-dependent optimizations for current CPU.")
if(CompileNative)
	#--- GCC architecture-dependent optimization
	check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
	if(HAS_MARCH_NATIVE)
		set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -march=native")
	endif()
	#--- Intel architecture-dependent optimization
	check_cxx_compiler_flag(-fast HAS_FAST)
	if(HAS_FAST)
		set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -fast")
	endif()
endif()

#MPI I/O options:
option(MPISafeWrite "Safe mode for MPI writes: write from one process; may be needed on some NFS installations")
if(MPISafeWrite)
	add_definitions("-DMPI_SAFE_WRITE")
endif()

#Workarounds for Windows compilation:
if(CYGWIN)
	add_definitions("-D_GNU_SOURCE")
	file(WRITE ${CMAKE_BINARY_DIR}/execinfo.h
"#ifndef EXECINFO_H
#define EXECINFO_H
inline int backtrace(void **buffer, int size) { return 0; }
inline char **backtrace_symbols(void *const *buffer, int size) { return 0; }
#endif //EXECINFO_H
"
	)
endif()

#------------- Macro to set JDFTx flags -------------

macro(set_JDFTx_flags targetName isGpu)
	if(${isGpu})
		set(GPU_ENABLED "-DGPU_ENABLED")
	else()
		set(GPU_ENABLED "")
	endif()
	set_target_properties(${targetName} PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS} ${GPU_ENABLED}"
		LINK_FLAGS "${EXTRA_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}"
		INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
		INSTALL_RPATH_USE_LINK_PATH TRUE)
endmacro()

#----------------------- Regular CPU targets ----------------

#External libraries to link to
set(EXTERNAL_LIBS ${HDF5_LIBRARIES} ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY} ${CBLAS_LAPACK_FFT_LIBRARIES} ${LIBXC_LIBRARY} ${EXTRA_LIBRARIES})

#Link options:
if(StaticLinking)
	set(LINK_TYPE STATIC)
	set(LINK_PREFIX -Wl,--whole-archive)
	set(LINK_SUFFIX -Wl,--no-whole-archive)
else()
	set(LINK_TYPE SHARED)
endif()

#Library with all the functionality:
FILE(GLOB jdftxlibSources core/*.cpp fluid/*.cpp electronic/*.cpp commands/*.cpp)
add_library(jdftxlib ${LINK_TYPE} ${jdftxlibSources})
target_link_libraries(jdftxlib ${EXTERNAL_LIBS})
set_target_properties(jdftxlib PROPERTIES OUTPUT_NAME "jdftx")
install(TARGETS jdftxlib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
set_JDFTx_flags(jdftxlib OFF)

#--------------- Nvidia GPU support using CUDA -----------------------
#----- Add another set of libraries and executables with a _gpu suffix
option(EnableCUDA "Enable NVIDIA GPU support using CUDA" OFF)
option(EnableCuSolver "Whether to use cuSolver GPU LAPACK (Requires CUDA >= 9)" OFF)

if(EnableCUDA)
	find_package(CUDA REQUIRED)
	set(CUDA_ARCH "compute_35" CACHE STRING "CUDA virtual architecture to compile for")
	set(CUDA_CODE "sm_35" CACHE STRING "CUDA gpu feature set version (sm_*) to compile for")
	set(CUDA_AUX_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})
	
	#remove libcuda.so from CUDA_LIBRARIES and save to CUDART_LIBRARY
	foreach(lib IN LISTS CUDA_LIBRARIES)
		if("${lib}" MATCHES "libcuda.so")
			list(REMOVE_ITEM CUDA_LIBRARIES "${lib}")
		endif("${lib}" MATCHES "libcuda.so")
	endforeach(lib)
	include_directories(${CUDA_INCLUDE_DIRS})

	#Check for cuSolver:
	if(EnableCuSolver)
		if(CUDA_VERSION_MAJOR LESS 9)
			message(FATAL_ERROR "cuSolver support requires CUDA >= 9.")
		endif()
		find_library(CUSOLVER_LIBRARY NAMES cusolver PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
		if(CUSOLVER_LIBRARY)
			find_package(OpenMP REQUIRED)
			set(CUDA_AUX_LIBRARIES ${CUDA_AUX_LIBRARIES} ${CUSOLVER_LIBRARY} ${OpenMP_CXX_FLAGS})
			add_definitions("-DCUSOLVER_ENABLED")
		else()
			message(FATAL_ERROR "Could not find the cuSolver library in the CUDA toolkkt located at ${CUDA_TOOLKIT_ROOT_DIR}.")
		endif()
	endif()
	message(STATUS "CUDA_LIBRARIES = ${CUDA_LIBRARIES};${CUDA_AUX_LIBRARIES}")
	
	set(CUDA_NVCC_FLAGS "-D_FORCE_INLINES;${CUDA_NVCC_FLAGS};-arch=${CUDA_ARCH};-code=${CUDA_CODE};-DGPU_ENABLED;--compiler-options;-fpic,-Wno-unused-but-set-variable")
	message(STATUS "CUDA_NVCC_FLAGS = ${CUDA_NVCC_FLAGS}")

	FILE(GLOB gpukernelsSources core/*.cu fluid/*.cu electronic/*.cu commands/*.cu tests/*.cu)
	cuda_add_library(gpukernels ${gpukernelsSources})
	target_link_libraries(gpukernels ${CUDA_AUX_LIBRARIES})

	#Library with all the functionality:
	FILE(GLOB jdftxlibSources core/*.cpp fluid/*.cpp electronic/*.cpp commands/*.cpp)
	add_library(jdftxlib_gpu ${LINK_TYPE} ${jdftxlibSources})
	target_link_libraries(jdftxlib_gpu ${EXTERNAL_LIBS} gpukernels)
	set_target_properties(jdftxlib_gpu PROPERTIES OUTPUT_NAME "jdftx_gpu")
	install(TARGETS jdftxlib_gpu LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
	set_JDFTx_flags(jdftxlib_gpu ON)
endif()

option(PinnedHostMemory "Use page-locked host memory to speed up GPU transfers" OFF)
if(PinnedHostMemory)
	add_definitions("-DPINNED_HOST_MEMORY")
endif()

option(CudaAwareMPI "Whether MPI library supports directly accessing CUDA pointer" OFF)
if(CudaAwareMPI)
	add_definitions("-DCUDA_AWARE_MPI")
endif()

#------------- Macro to add an executable pair (regular and GPU) -------------

macro(add_JDFTx_executable execName execSources)
	#CPU version
	add_executable(${execName} ${ARGN} ${execSources})
	target_link_libraries(${execName} ${LINK_PREFIX} jdftxlib ${LINK_SUFFIX})
	if(NOT "${ARGN}" MATCHES EXCLUDE_FROM_ALL)
		install(TARGETS ${execName} RUNTIME DESTINATION bin)
	endif()
	set_JDFTX_flags(${execName} OFF)
	
	if(EnableCUDA) #GPU version
		add_executable(${execName}_gpu ${ARGN} ${execSources})
		target_link_libraries(${execName}_gpu ${LINK_PREFIX} jdftxlib_gpu ${LINK_SUFFIX})
		if(NOT "${ARGN}" MATCHES EXCLUDE_FROM_ALL)
			install(TARGETS ${execName}_gpu RUNTIME DESTINATION bin)
		endif()
		set_JDFTx_flags(${execName}_gpu ON)
	endif()
endmacro()

#----------------------- Primary interface executable ------------------------

add_JDFTx_executable(jdftx jdftx.cpp)

#----------------------- Secondary executables ------------------------

#Wannier:
FILE(GLOB wannierSources wannier/*.cpp)
add_JDFTx_executable(wannier "${wannierSources}")

#Phonon:
FILE(GLOB phononSources phonon/*.cpp)
add_JDFTx_executable(phonon "${phononSources}")

#-----------------------------------------------------------------------------

#Documentation via Doxygen:
find_package(Doxygen)
if(DOXYGEN_FOUND)
	set(DOC_TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/doc.tmp")
	file(MAKE_DIRECTORY ${DOC_TMP_DIR})
	file(MAKE_DIRECTORY ${DOC_TMP_DIR}/commands)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)
	#Target for running doxygen:
	add_custom_target(doc
		COMMAND echo "--- Autogenerating command documentation"
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/jdftx   -w jdftx   > "${DOC_TMP_DIR}/commands/jdftx.dox"
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/wannier -w wannier > "${DOC_TMP_DIR}/commands/wannier.dox"
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/phonon  -w phonon  > "${DOC_TMP_DIR}/commands/phonon.dox"
		COMMAND echo "--- Updating command index"
		COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/doc/commands/category_update.sh "${DOC_TMP_DIR}" > "${DOC_TMP_DIR}/commands_categorized.dox"
		COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/doc/commands/index_update.sh "${DOC_TMP_DIR}" > "${DOC_TMP_DIR}/commands_index.dox"
		COMMAND echo "--- Updating scripts documentation"
		COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/doc/scripts/index_update.sh "${CMAKE_CURRENT_SOURCE_DIR}" > "${DOC_TMP_DIR}/scripts_manual.dox"
		COMMAND echo "--- Running doxygen"
		COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/doc/html
		COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		DEPENDS jdftx wannier phonon
	)
endif()

#Tests
enable_testing()
add_subdirectory(test)

#Optional: additional auxiliary test executables, mostly related to fluid development
add_subdirectory(aux)

#Optional features: primarily interfaces to other codes
add_subdirectory(opt)

#Pseudopotential library
add_subdirectory(pseudopotentials)


#----------------------- Installation --------------------------
#Executables and library install rules already added above
#--- headers
install(DIRECTORY
	commands
	core
	electronic
	fluid
	phonon
	wannier
	DESTINATION include/jdftx
	FILES_MATCHING PATTERN "*.h")
install(FILES "${CMAKE_BINARY_DIR}/config.h" DESTINATION include/jdftx)
#--- pseudopotentials
install(DIRECTORY
	"${CMAKE_BINARY_DIR}/pseudopotentials/GBRV"
	"${CMAKE_BINARY_DIR}/pseudopotentials/SG15"
	DESTINATION share/jdftx/pseudopotentials)
#--- scripts
install(DIRECTORY scripts DESTINATION share/jdftx USE_SOURCE_PERMISSIONS)
