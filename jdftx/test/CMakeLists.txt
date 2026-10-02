add_custom_target(testresults COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/printResults.sh ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_target(testclean COMMAND rm -f */*.out */*.wfns */*.fillings */*.ionpos */*.eigenvals */*.fluidState */results */summary WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

macro(add_jdftx_test testName)
	add_test(NAME ${testName} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh ${testName} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR})
endmacro()

add_jdftx_test(openShell)
add_jdftx_test(vibrations)
add_jdftx_test(moleculeSolvation)
add_jdftx_test(ionSolvation)
add_jdftx_test(latticeOpt)
add_jdftx_test(metalBulk)
add_jdftx_test(plusU)
add_jdftx_test(spinOrbit)
add_jdftx_test(graphene)
add_jdftx_test(metalSurface)
add_jdftx_test(checkpoint)

add_jdftx_test(defectCoulomb)

# End-to-end one-shot generation, native potential comparison and an NSCF solve.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME defectPotential COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/defectCoulomb/checkPotential.py
        ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/defectCoulomb
        ${CMAKE_SOURCE_DIR}/scripts/defectDensityToVscloc)
    set_tests_properties(defectCoulomb PROPERTIES FIXTURES_SETUP defectDensities)
    set_tests_properties(defectPotential PROPERTIES FIXTURES_REQUIRED defectDensities)
    add_test(NAME defectEmbedding COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/defectCoulomb/checkEmbedding.py
        ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/defectCoulomb
        ${CMAKE_SOURCE_DIR}/scripts/defectDensityToVscloc)
endif()
