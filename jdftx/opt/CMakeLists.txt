option(EnableVASP "If yes, create a shared library exposing solvation models to VASP.")

if(EnableVASP)
	add_subdirectory(VASPinterface)
endif()