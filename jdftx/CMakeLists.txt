project(dft)

cmake_minimum_required(VERSION 2.6)
enable_language (Fortran)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake-Modules/")

#Package configuration:
set(CPACK_PACKAGE_NAME "JDFTx")
set(CPACK_PACKAGE_DESCRIPTION "The playground for joint density functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "R&D code from the Arias group @ Cornell Physics")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "preAlpha")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

#Get the subversion revision number
find_package(Subversion)
if(Subversion_FOUND)
        Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
        set(SUBVERSION_REVISION ${ER_WC_REVISION})
        message(STATUS "SVN revision number: ${SUBVERSION_REVISION}")
endif(Subversion_FOUND)

find_package(LAPACK REQUIRED)

find_package(CBLAS REQUIRED)
include_directories(${CBLAS_INCLUDE_DIR})

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

find_package(LIBXC)
if(LIBXC_FOUND)
	include_directories(${LIBXC_INCLUDE_DIR})
	add_definitions("-DLIBXC_ENABLED")
else(LIBXC_FOUND)
	set(LIBXC_LIBRARY "")
endif(LIBXC_FOUND)

#Process configuration information into config.h (with config.in.h as a template)
configure_file(${CMAKE_SOURCE_DIR}/config.in.h ${CMAKE_BINARY_DIR}/config.h)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR})

option(EnableProfiling "Optionally enable profiling to determine ideal functions to optimize (for developers)")
if(EnableProfiling)
	add_definitions("-DENABLE_PROFILING")
endif(EnableProfiling)

#Compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-result -g -O3")

#----------------------- Regular CPU stuff ----------------

#Set of libraries used by each target
set(CORE_LIBS ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${CBLAS_LIBRARY} ${GSL_LIBRARY} ${LAPACK_LIBRARIES} ${LIBXC_LIBRARY})

#Core Library:
FILE(GLOB coreSources core/*.cpp)
add_library(jdftxCore SHARED ${coreSources})
target_link_libraries(jdftxCore ${CORE_LIBS})

#Fluid library:
FILE(GLOB fluidSources fluid/*.cpp)
add_library(jdftxFluid SHARED ${fluidSources})
target_link_libraries(jdftxFluid jdftxCore)

#Electronic library:
FILE(GLOB electronicSources electronic/*.cpp)
add_library(jdftxElectronic SHARED ${electronicSources})
target_link_libraries(jdftxElectronic jdftxFluid)

#Commands and parser:
FILE(GLOB commandSources commands/*.cpp)
add_library(jdftxCommands SHARED ${commandSources})
target_link_libraries(jdftxCommands jdftxElectronic)

#Main executable:
add_executable(jdftx jdftx.cpp)
target_link_libraries(jdftx jdftxCommands)

#Target properties:
set_target_properties(jdftx jdftxCore jdftxFluid jdftxElectronic jdftxCommands  PROPERTIES COMPILE_FLAGS "-std=c++0x")

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

	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-arch=sm_20;-DGPU_ENABLED") #Compile for Fermi or higher (change this for older cards)

	FILE(GLOB gpukernelsSources core/*.cu fluid/*.cu electronic/*.cu commands/*.cu tests/*.cu)
	cuda_add_library(gpukernels SHARED ${gpukernelsSources})
	target_link_libraries(gpukernels ${CUDA_CUBLAS_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})

	#Core Library (Gpu Version):
	FILE(GLOB coreSources core/*.cpp)
	add_library(jdftxCore_gpu SHARED ${coreSources})
	target_link_libraries(jdftxCore_gpu ${CORE_LIBS} gpukernels)

	#Fluid library (Gpu Version):
	FILE(GLOB fluidSources fluid/*.cpp)
	add_library(jdftxFluid_gpu SHARED ${fluidSources})
	target_link_libraries(jdftxFluid_gpu jdftxCore_gpu)

	#Electronic library (Gpu Version):
	FILE(GLOB electronicSources electronic/*.cpp)
	add_library(jdftxElectronic_gpu SHARED ${electronicSources})
	target_link_libraries(jdftxElectronic_gpu jdftxFluid_gpu)

	#Commands and parser (Gpu version):
	FILE(GLOB commandSources commands/*.cpp)
	add_library(jdftxCommands_gpu SHARED ${commandSources})
	target_link_libraries(jdftxCommands_gpu jdftxElectronic_gpu)

	#Main executable (Gpu version):
	add_executable(jdftx_gpu jdftx.cpp)
	target_link_libraries(jdftx_gpu jdftxCommands_gpu)

	#Target properties for Gpu versions:
	set_target_properties(jdftx_gpu jdftxCore_gpu jdftxFluid_gpu jdftxElectronic_gpu jdftxCommands_gpu PROPERTIES COMPILE_FLAGS "-std=c++0x -DGPU_ENABLED")

endif(CUDA_FOUND)

#-----------------------------------------------------------------------------


#Documentation via Doxygen:
find_package(Doxygen)
if(DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        add_custom_target(doc
                ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif (DOXYGEN_FOUND)

#Additional test executables, mostly related to the fluid
add_subdirectory(tests)
