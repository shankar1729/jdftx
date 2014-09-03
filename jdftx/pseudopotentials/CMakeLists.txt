#Pseudopotential library target
add_custom_target(PseudopotentialLibrary ALL
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/update.sh ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
