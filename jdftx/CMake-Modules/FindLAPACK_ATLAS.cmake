#Wrapper to system FindLAPACK.cmake that handles new ATLAS library format correctly

if(LAPACK_LIBRARIES)
	#Externally specified LAPACK library list (with dependencies such as BLAS)
	message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
else()
	#Try to find the ATLAS sequential library:
	find_library(ATLAS_LIBRARY NAMES satlas PATHS ${ATLAS_PATH} ${ATLAS_PATH}/lib ${ATLAS_PATH}/lib64 /usr/lib64/atlas /usr/lib/atlas NO_DEFAULT_PATH)
	find_library(ATLAS_LIBRARY NAMES satlas)

	if(ATLAS_LIBRARY)
		set(LAPACK_LIBRARIES ${ATLAS_LIBRARY})
		set(LAPACK_FOUND TRUE)
		if(NOT LAPACK_FIND_QUIETLY)
			message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
		endif()
	else()
		#If not found, invoke the system FindLAPACK:
		if(LAPACK_ATLAS_FIND_REQUIRED)
			find_package(LAPACK REQUIRED)
		else()
			find_package(LAPACK)
		endif()
	endif()
endif()

if(EnableScaLAPACK)
	#Try to find the ScaLAPACK library
	find_library(SCALAPACK_LIBRARY NAMES scalapack PATHS ${SCALAPACK_PATH} ${SCALAPACK_PATH}/lib ${SCALAPACK_PATH}/lib64 NO_DEFAULT_PATH)
	find_library(SCALAPACK_LIBRARY NAMES scalapack)
	if(SCALAPACK_LIBRARY)
		message(STATUS "Found ScaLAPACK: ${SCALAPACK_LIBRARY}")
		set(LAPACK_LIBRARIES ${SCALAPACK_LIBRARY} ${LAPACK_LIBRARIES})
		add_definitions("-DSCALAPACK_ENABLED")
	else()
		message(FATAL_ERROR "Could not find the ScaLAPACK library (Add -D SCALAPACK_PATH=<path> to the cmake commandline for a non-standard installation)")
	endif()
endif()
