add_library(jdftxVASPinterface SHARED FluidInterface.cpp)
target_link_libraries(jdftxVASPinterface jdftxlib)

add_executable(test test.cpp)
target_link_libraries(test jdftxVASPinterface)

set_target_properties(jdftxVASPinterface test
	PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS}"
		LINK_FLAGS "${EXTRA_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}")

