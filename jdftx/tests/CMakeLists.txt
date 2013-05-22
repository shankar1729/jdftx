# Various test and auxiliary utility codes

set(targetNameList
	TestRandom          #Random tests during develoment (no specific purpose)
	TestSphere          #Hard sphere (surface tension)
	Capacitor           #Parallel plate capacitor (dielectric constant)
	TestPlanar          #Liquid-Vapor interface (surface tension)
	TestOperators       #Test operators and memory management
	SO3quadConvergence  #Test SO3 quadrature convergence
	NonlinearEps        #Nonlinear dielectric constant
	CavitationH2O       #Compute cavitation energy given electron density (to complement JDFT1)
	CavitationLJ        #Compute hard sphere cavitation energy given electron density
	TestFixedN          #Testing fixed N (canonical) classical DFT capability
	TestGaussian        #Tests water functionals with parabolic potential well
	SolvatedCMOS        #Floating gate transistor in ionic liquid
	SphericalChi        #Compute spherical decomposition of non-local susceptibility
	ElectrostaticRadius #Estimate electrostatic radius of solvent molecule
)

foreach(targetName ${targetNameList})
	#CPU version
	add_executable(${targetName} ${targetName}.cpp)
	target_link_libraries(${targetName} jdftxlib)
	set_target_properties(${targetName} PROPERTIES
		COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS}"
		LINK_FLAGS "${EXTRA_CXX_FLAGS}")
	
	if(CUDA_FOUND) #GPU version
	        add_executable(${targetName}_gpu ${targetName}.cpp)
	        target_link_libraries(${targetName}_gpu jdftxlib_gpu)
	        set_target_properties(${targetName}_gpu PROPERTIES
	                COMPILE_FLAGS "${EXTRA_CXX_FLAGS} ${JDFTX_CPU_FLAGS} -DGPU_ENABLED"
	                LINK_FLAGS "${EXTRA_CXX_FLAGS}")
	endif()
endforeach()

