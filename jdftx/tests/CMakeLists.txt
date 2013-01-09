# Various test codes, mostly from the original Water3D code

ADD_EXECUTABLE(         TestRandom            TestRandom.cpp) #Random tests during develoment (no specific purpose)
ADD_EXECUTABLE(         TestSphere            TestSphere.cpp) #Hard sphere (surface tension)
ADD_EXECUTABLE(         TestPlanar            TestPlanar.cpp) #Parallel plate capacitor (dielectric constant)
ADD_EXECUTABLE(      TestOperators         TestOperators.cpp) #Test operators and memory management
ADD_EXECUTABLE( SO3quadConvergence    SO3quadConvergence.cpp) #Test SO3 quadrature convergence
ADD_EXECUTABLE(       NonlinearEps          NonlinearEps.cpp) #Nonlinear dielectric constant
ADD_EXECUTABLE(      CavitationH2O         CavitationH2O.cpp) #Compute cavitation energy given electron density (to complement JDFT1)
ADD_EXECUTABLE(    CavitationLJ             CavitationLJ.cpp) #Compute hard sphere cavitation energy given electron density
ADD_EXECUTABLE(         TestFixedN            TestFixedN.cpp) #Attemping Fixed N capability
ADD_EXECUTABLE(       TestGaussian          TestGaussian.cpp) #Tests water functionals with parabolic potential well
ADD_EXECUTABLE(       SolvatedCMOS          SolvatedCMOS.cpp) #Floating gate transistor in ionic liquid
TARGET_LINK_LIBRARIES(         TestRandom   jdftxFluid)
TARGET_LINK_LIBRARIES(         TestSphere   jdftxFluid)
TARGET_LINK_LIBRARIES(         TestPlanar   jdftxFluid)
TARGET_LINK_LIBRARIES(      TestOperators   jdftxElectronic)
TARGET_LINK_LIBRARIES( SO3quadConvergence   jdftxFluid)
TARGET_LINK_LIBRARIES(       NonlinearEps   jdftxFluid)
TARGET_LINK_LIBRARIES(      CavitationH2O   jdftxFluid)
TARGET_LINK_LIBRARIES(    CavitationLJ      jdftxFluid)
TARGET_LINK_LIBRARIES(         TestFixedN   jdftxFluid)
TARGET_LINK_LIBRARIES(       TestGaussian   jdftxFluid)
TARGET_LINK_LIBRARIES(       SolvatedCMOS   jdftxFluid)

set_target_properties(
               TestRandom
               TestSphere
               TestPlanar
            TestOperators
       SO3quadConvergence
             NonlinearEps
            CavitationH2O
             CavitationLJ
               TestFixedN
             TestGaussian
             SolvatedCMOS
	PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS}"
		LINK_FLAGS "${EXTRA_CXX_FLAGS}")

#---------------------------------------------------------------------------------------------------
if(CUDA_FOUND)
	#Copies of the above, but now with the GPU flag enabled
	#All kernels would have been compiled into the gpukernels target in the parent CMakeLists.txt
	#These are regular C++0x (not cuda) builds, just with GPU_ENABLED defined

	ADD_EXECUTABLE(         TestRandom_gpu            TestRandom.cpp) #Random tests during develoment (no specific purpose)
	ADD_EXECUTABLE(         TestSphere_gpu            TestSphere.cpp) #Hard sphere (surface tension)
	ADD_EXECUTABLE(         TestPlanar_gpu            TestPlanar.cpp) #Parallel plate capacitor (dielectric constant)
	ADD_EXECUTABLE(      TestOperators_gpu         TestOperators.cpp) #Test operators and memory management
	ADD_EXECUTABLE( SO3quadConvergence_gpu    SO3quadConvergence.cpp) #Test operators and memory management
	ADD_EXECUTABLE(       NonlinearEps_gpu          NonlinearEps.cpp) #Nonlinear dielectric constant
	ADD_EXECUTABLE(      CavitationH2O_gpu         CavitationH2O.cpp) #Compute cavitation energy given electron density (to complement JDFT1)
	ADD_EXECUTABLE(       CavitationLJ_gpu          CavitationLJ.cpp) #Compute cavitation energy given electron density
	ADD_EXECUTABLE(         TestFixedN_gpu            TestFixedN.cpp) #Attemping Fixed N capability
	ADD_EXECUTABLE(       TestGaussian_gpu          TestGaussian.cpp) #Tests water functionals with parabolic potential well
	ADD_EXECUTABLE(       SolvatedCMOS_gpu          SolvatedCMOS.cpp) #Floating gate transistor in ionic liquid

	TARGET_LINK_LIBRARIES(         TestRandom_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(         TestSphere_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(         TestPlanar_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(      TestOperators_gpu   jdftxElectronic_gpu)
	TARGET_LINK_LIBRARIES( SO3quadConvergence_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(       NonlinearEps_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(      CavitationH2O_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(       CavitationLJ_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(         TestFixedN_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(       TestGaussian_gpu   jdftxFluid_gpu)
	TARGET_LINK_LIBRARIES(       SolvatedCMOS_gpu   jdftxFluid_gpu)

	set_target_properties(
	               TestRandom_gpu
	               TestSphere_gpu
	               TestPlanar_gpu
	            TestOperators_gpu
	       SO3quadConvergence_gpu
	             NonlinearEps_gpu
	            CavitationH2O_gpu
                 CavitationLJ_gpu
	               TestFixedN_gpu
	             TestGaussian_gpu
                 SolvatedCMOS_gpu
	PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS} -DGPU_ENABLED"
		LINK_FLAGS "${EXTRA_CXX_FLAGS}")

endif()

