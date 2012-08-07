project(fluid1D)
cmake_minimum_required(VERSION 2.6)

set(CPACK_PACKAGE_NAME "fluid1D")
set(CPACK_PACKAGE_DESCRIPTION "One-dimensional classical density functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "R&D code from the Arias group @ Cornell Physics")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "preAlpha")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

add_definitions(-Wall -g -O3 -std=c++0x)

set(JDFTX_SRC "<full-path-to-jdftx>" CACHE STRING "Full path to the jdftx source root")
set(CMAKE_MODULE_PATH "${JDFTX_SRC}/CMake-Modules/")

#Get the subversion revision number
find_package(Subversion)
if(Subversion_FOUND)
        Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
        set(SUBVERSION_REVISION ${ER_WC_REVISION})
        message(STATUS "SVN revision number: ${SUBVERSION_REVISION}")
endif(Subversion_FOUND)

find_package(CBLAS REQUIRED)
include_directories(${CBLAS_INCLUDE_DIR})

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

include_directories("${JDFTX_SRC}")

#Process configuration information into config.h (with config.in.h as a template)
configure_file(${CMAKE_SOURCE_DIR}/config.in.h ${CMAKE_BINARY_DIR}/config.h)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR})

#Core 1D library:
FILE(GLOB coreSources core1D/*.cpp)
add_library(core1D SHARED ${coreSources}
	${JDFTX_SRC}/core/Util.cpp ${JDFTX_SRC}/core/AutoThreadCount.cpp
	${JDFTX_SRC}/core/BlasExtra.cpp ${JDFTX_SRC}/core/Random.cpp 
)
target_link_libraries(core1D ${FFTW3_THREADS_LIBRARY} ${FFTW3_LIBRARY} ${GSL_LIBRARY} ${CBLAS_LIBRARY})

#Fluid library:
FILE(GLOB fluidSources fluid1D/*.cpp)
add_library(fluid1D SHARED ${fluidSources})
target_link_libraries(fluid1D core1D)

#Executables
add_subdirectory(tests)
