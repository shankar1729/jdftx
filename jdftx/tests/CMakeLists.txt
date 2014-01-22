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
	SigmaModel          #Surface tension from classical DFT and weighted-density model
	SolvatedCMOS        #Floating gate transistor in ionic liquid
	SphericalChi        #Compute spherical decomposition of non-local susceptibility
	ElectrostaticRadius #Estimate electrostatic radius of solvent molecule
	SlaterDetOverlap    #Estimate the dipole matrix element of two column bundles
)

foreach(targetName ${targetNameList})
	add_JDFTx_executable(${targetName} ${targetName}.cpp)
endforeach()
