
add_custom_target(testresults COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/printResults.sh ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_target(testclean COMMAND rm -f */*.out */*.wfns */*.fillings */*.ionpos */*.Haux */*.fluidState */results */summary WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

macro(add_jdftx_test testName)
	add_test(NAME ${testName} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh ${testName} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR})
endmacro()

add_jdftx_test(vibrations)
add_jdftx_test(moleculeSolvation)