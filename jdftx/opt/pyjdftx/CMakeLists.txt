cmake_minimum_required(VERSION 3.15...4.0)
project(pyjdftx LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

execute_process(
	COMMAND python3 -c "import mpi4py; print(mpi4py.get_include())"
	OUTPUT_VARIABLE MPI4PY_INCLUDE
)
string(STRIP ${MPI4PY_INCLUDE} MPI4PY_INCLUDE)
message(STATUS "Found mpi4py includes: ${MPI4PY_INCLUDE}")

find_library(JDFTX_LIBRARY NAMES jdftx)
if(JDFTX_LIBRARY)
	message(STATUS "Found libjdftx: ${JDFTX_LIBRARY}")
else()
	message(FATAL_ERROR "Could not find jdftx library (check LIBRARY_PATH)")
endif()

set(INCLUDE_DIRS ${MPI4PY_INCLUDE} ${CMAKE_SOURCE_DIR}/../..)

pybind11_add_module(pyjdftx src/pyjdftx.cpp)
target_include_directories(pyjdftx PRIVATE ${INCLUDE_DIRS})
target_link_libraries(pyjdftx PRIVATE ${JDFTX_LIBRARY})
target_compile_definitions(pyjdftx PRIVATE MPI_ENABLED JDFTX_LIBRARY="\\\"${JDFTX_LIBRARY}\\\"")
install(TARGETS pyjdftx DESTINATION pyjdftx)

find_library(JDFTX_GPU_LIBRARY NAMES jdftx_gpu)
if(JDFTX_GPU_LIBRARY)
	message(STATUS "Found libjdftx_gpu: ${JDFTX_LIBRARY} (will build gpu version)")
	pybind11_add_module(pyjdftx_gpu src/pyjdftx.cpp)
	target_include_directories(pyjdftx_gpu PRIVATE ${INCLUDE_DIRS})
	target_link_libraries(pyjdftx_gpu PRIVATE ${JDFTX_GPU_LIBRARY})
	target_compile_definitions(pyjdftx_gpu PRIVATE MPI_ENABLED GPU_ENABLED JDFTX_LIBRARY="\\\"${JDFTX_GPU_LIBRARY}\\\"")
	install(TARGETS pyjdftx_gpu DESTINATION pyjdftx)
endif()
