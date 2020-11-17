# Various  auxiliary test and utility codes

set(targetNameList
	TestRandom          #Random tests during develoment (no specific purpose)
	TestSphere          #Hard sphere (surface tension)
	Capacitor           #Parallel plate capacitor (dielectric constant)
	TestPlanar          #Liquid-Vapor interface (surface tension)
	TestOperators       #Test operators and memory management
	SO3quadConvergence  #Test SO3 quadrature convergence
	NonlinearEps        #Nonlinear dielectric constant
	TestGaussian        #Tests water functionals with parabolic potential well
	SigmaModel          #Surface tension from classical DFT and weighted-density model
	SphericalChi        #Compute spherical decomposition of non-local susceptibility
	ElectrostaticRadius #Estimate electrostatic radius of solvent molecule
	SlaterDetOverlap    #Estimate the dipole matrix element of two column bundles
	TestSchrodinger     #Davidson solution of Schrodinger equation as a performance benchmark
)

foreach(targetName ${targetNameList})
	add_JDFTx_executable(${targetName} ${targetName}.cpp EXCLUDE_FROM_ALL)
endforeach()
add_custom_target(aux DEPENDS ${targetNameList})

