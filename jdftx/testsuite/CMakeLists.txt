
add_custom_target(testresults COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/printResults.sh ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_target(testclean COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*/*.out ${CMAKE_CURRENT_BINARY_DIR}/*/results ${CMAKE_CURRENT_BINARY_DIR}/*/summary )

macro(add_jdftx_test testName)
	add_test(NAME ${testName} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh ${testName} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR})
endmacro()

add_jdftx_test(vibrations)
