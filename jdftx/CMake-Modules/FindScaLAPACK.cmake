#Check for built in scalapack:
find_library(SCALAPACK_LIBRARY NAMES scalapack PATHS ${SCALAPACK_PATH} ${SCALAPACK_PATH}/lib ${SCALAPACK_PATH}/lib64 NO_DEFAULT_PATH)
find_library(SCALAPACK_LIBRARY NAMES scalapack)

if(SCALAPACK_LIBRARY)
	set(SCALAPACK_FOUND TRUE)
endif()

if(SCALAPACK_FOUND)
	if(NOT SCALAPACK_FIND_QUIETLY)
		message(STATUS "Found ScaLAPACK: ${SCALAPACK_LIBRARY}")
	endif()
else()
	if(ScaLAPACK_FIND_REQUIRED)
		if(NOT SCALAPACK_LIBRARY)
			message(FATAL_ERROR "Could not find the SCALAPACK shared library. (Add -D SCALAPACK_PATH=<path> to the cmake commandline for a non-standard installation)")
		endif()
	endif()
endif()
