project(dft)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake-Modules/")

#Package configuration:
set(CPACK_PACKAGE_NAME "JDFTx")
set(CPACK_PACKAGE_DESCRIPTION "The playground for joint density functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "R&D code from the Arias group @ Cornell Physics")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "99")
set(CPACK_PACKAGE_VERSION_PATCH "alpha")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

#Get the subversion revision number
find_package(Subversion)
if(Subversion_FOUND)
        Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
        set(SUBVERSION_REVISION ${ER_WC_REVISION})
        message(STATUS "SVN revision number: ${SUBVERSION_REVISION}")
endif()

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

#Use indepedent CBLAS, LAPACK and FFTW (default) or MKL for all three
option(EnableMKL "Use Intel MKL to provide BLAS, LAPACK and FFTs")
option(ForceFFTW "Force usage of FFTW (even if MKL is enabled)")
if(EnableMKL)
	find_package(MKL REQUIRED)
	set(CBLAS_LAPACK_FFT_LIBRARIES ${MKL_LIBRARIES}) #MKL provides CBLAS, FFTW3 and LAPACK
	add_definitions("-DMKL_PROVIDES_BLAS") #Special handling is required for BLAS calls from threaded sections
	include_directories(${MKL_INCLUDE_DIR})
	if(ForceFFTW)
		find_package(FFTW3 REQUIRED)
		set(CBLAS_LAPACK_FFT_LIBRARIES ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${MKL_LIBRARIES}) #Explicit FFTW3, rest from MKL
	else()
		add_definitions("-DMKL_PROVIDES_FFT") #Special handling is required for FFT initialization
		set(CBLAS_LAPACK_FFT_LIBRARIES ${MKL_LIBRARIES}) #MKL provides CBLAS, FFTW3 and LAPACK
	endif()
else()
	find_package(FFTW3 REQUIRED)
	find_package(LAPACK REQUIRED)
	find_package(CBLAS REQUIRED)
	set(CBLAS_LAPACK_FFT_LIBRARIES ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${CBLAS_LIBRARY} ${LAPACK_LIBRARIES})
endif()
include_directories(${FFTW3_INCLUDE_DIR})

option(EnableMPI "Use MPI parallelization (in addition to threads / gpu)" ON)
if(EnableMPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	add_definitions("-DMPI_ENABLED")
endif()

option(EnableLibXC "Use LibXC to provide additional exchange-correlation functionals")
if(EnableLibXC)
	find_package(LIBXC REQUIRED)
	include_directories(${LIBXC_INCLUDE_DIR})
	add_definitions("-DLIBXC_ENABLED")
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
#--- Link-time optimizations (GCC)
check_cxx_compiler_flag(-flto HAS_FLTO)
if(HAS_FLTO)
	set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -flto")
endif()
#--- Link-time optimizations (ICC)
check_cxx_compiler_flag(-ipo HAS_IPO)
if(HAS_IPO)
	set(JDFTX_CPU_FLAGS "${JDFTX_CPU_FLAGS} -ipo")
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

#------------- Macro to set JDFTx flags -------------

macro(set_JDFTx_flags targetName isGpu)
	if(${isGpu})
		set(GPU_ENABLED "-DGPU_ENABLED")
	else()
		set(GPU_ENABLED "")
	endif()
	set_target_properties(${targetName} PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS} ${GPU_ENABLED}"
		LINK_FLAGS "${EXTRA_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endmacro()

#----------------------- Regular CPU targets ----------------

#External libraries to link to
set(EXTERNAL_LIBS ${MPI_CXX_LIBRARIES} ${CBLAS_LAPACK_FFT_LIBRARIES} ${GSL_LIBRARY} ${LIBXC_LIBRARY})

#Library with all the functionality:
FILE(GLOB jdftxlibSources core/*.cpp fluid/*.cpp electronic/*.cpp commands/*.cpp)
add_library(jdftxlib SHARED ${jdftxlibSources})
target_link_libraries(jdftxlib ${EXTERNAL_LIBS})
set_target_properties(jdftxlib PROPERTIES OUTPUT_NAME "jdftx")
set_JDFTx_flags(jdftxlib OFF)

#--------------- Nvidia GPU support using CUDA -----------------------
#----- Add another set of libraries and executables with a _gpu suffix
find_package(CUDA)
if(CUDA_FOUND)
	#remove libcuda.so from CUDA_LIBRARIES and save to CUDART_LIBRARY
	foreach(lib IN LISTS CUDA_LIBRARIES)
		if("${lib}" MATCHES "libcuda.so")
			list(REMOVE_ITEM CUDA_LIBRARIES "${lib}")
		endif("${lib}" MATCHES "libcuda.so")
	endforeach(lib)
	message(STATUS "CUDA_LIBRARIES = ${CUDA_LIBRARIES}")
	include_directories(${CUDA_INCLUDE_DIRS})

	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-arch=compute_20;-code=sm_20,sm_30;-DGPU_ENABLED;--compiler-options;-fpic,-Wno-unused-but-set-variable") #Compile for Fermi or higher (change this for older cards)

	FILE(GLOB gpukernelsSources core/*.cu fluid/*.cu electronic/*.cu commands/*.cu tests/*.cu)
	cuda_add_library(gpukernels ${gpukernelsSources})
	target_link_libraries(gpukernels ${CUDA_CUBLAS_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})

	#Library with all the functionality:
	FILE(GLOB jdftxlibSources core/*.cpp fluid/*.cpp electronic/*.cpp commands/*.cpp)
	add_library(jdftxlib_gpu SHARED ${jdftxlibSources})
	target_link_libraries(jdftxlib_gpu ${EXTERNAL_LIBS} gpukernels)
	set_target_properties(jdftxlib_gpu PROPERTIES OUTPUT_NAME "jdftx_gpu")
	set_JDFTx_flags(jdftxlib_gpu ON)
endif()

#------------- Macro to add an executable pair (regular and GPU) -------------

macro(add_JDFTx_executable execName execSources)
	#CPU version
	add_executable(${execName} ${execSources})
	target_link_libraries(${execName} jdftxlib)
	set_JDFTX_flags(${execName} OFF)
	
	if(CUDA_FOUND) #GPU version
	        add_executable(${execName}_gpu ${execSources})
	        target_link_libraries(${execName}_gpu jdftxlib_gpu)
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
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)
	#Target for running doxygen:
	add_custom_target(doc
		COMMAND echo "--- Autogenerating command documentation"
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/jdftx   -w jdftx   > ${CMAKE_CURRENT_SOURCE_DIR}/commands/manual.dox
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/wannier -w wannier > ${CMAKE_CURRENT_SOURCE_DIR}/wannier/manual.dox
		COMMAND ${CMAKE_CURRENT_BINARY_DIR}/phonon  -w phonon  > ${CMAKE_CURRENT_SOURCE_DIR}/phonon/manual.dox
		COMMAND echo "--- Updating command index"
		COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/doc/commands/index_update.sh "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/doc/commands/index.dox"
		COMMAND echo "--- Running doxygen"
		COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		DEPENDS jdftx wannier phonon
	)
endif()

#Additional test executables, mostly related to the fluid
add_subdirectory(tests)

#Test suite
enable_testing()
add_subdirectory(testsuite)

#Optional features: primarily interfaces to other codes
add_subdirectory(opt)

#Pseudopotential library
add_subdirectory(pseudopotentials)
