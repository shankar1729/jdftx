#Wrapper to system FindLAPACK.cmake that handles new ATLAS library format correctly

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
	if(LAPACK_ATLAS_REQUIRED)
		find_package(LAPACK REQUIRED)
	else()
		find_package(LAPACK)
	endif()
endif()
