#include <electronic/VanDerWaalsD3.h>


VanDerWaalsD3::VanDerWaalsD3(const Everything& e) : VanDerWaals(e)
{
	logPrintf("\nInitializing DFT-D3 calculator:\n");

	Citations::add("DFT-D3 dispersion correction", "S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132, 154104 (2010)");
}


VanDerWaalsD3::~VanDerWaalsD3()
{
}


double VanDerWaalsD3::getScaleFactor(string exCorrName, double scaleOverride) const
{	return 0.;  //global scale factor not used in D3
}


double VanDerWaalsD3::energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT) const
{	die("Not yet implemented.\n");
	return 0.;
}
