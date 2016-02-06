project(fluid1D)
cmake_minimum_required(VERSION 2.6)

set(CPACK_PACKAGE_NAME "fluid1D")
set(CPACK_PACKAGE_DESCRIPTION "One-dimensional classical density functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Companion code to Comp. Phys. Comm. 185, 818 (2014)")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

add_definitions(-Wall -g -O3 -std=c++0x)

set(JDFTX_SRC "<full-path-to-jdftx>" CACHE STRING "Full path to the jdftx source root")
set(CMAKE_MODULE_PATH "${JDFTX_SRC}/CMake-Modules/")

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
add_library(fluid1D SHARED ${fluidSources}
	${JDFTX_SRC}/fluid/SO3quad.cpp ${JDFTX_SRC}/fluid/S2quad.cpp
)
target_link_libraries(fluid1D core1D)

#Executables
add_subdirectory(tests)
