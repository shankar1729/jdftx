add_executable(TestOperators TestOperators.cpp)
target_link_libraries(TestOperators core1D)

add_executable(TestSphere TestSphere.cpp)
target_link_libraries(TestSphere fluid1D)

add_executable(TestPlanar TestPlanar.cpp)
target_link_libraries(TestPlanar fluid1D)

add_executable(SigmaVsRadius SigmaVsRadius.cpp)
target_link_libraries(SigmaVsRadius fluid1D)

add_executable(TestRandom TestRandom.cpp)
target_link_libraries(TestRandom fluid1D)

add_executable(NonlinearEps NonlinearEps.cpp)
target_link_libraries(NonlinearEps fluid1D)

add_executable(CorrFunc CorrFunc.cpp)
target_link_libraries(CorrFunc fluid1D)

add_executable(TestCFC TestCFC.cpp)
target_link_libraries(TestCFC fluid1D)
