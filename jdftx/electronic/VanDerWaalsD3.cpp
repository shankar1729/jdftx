#include <electronic/VanDerWaalsD3.h>
#include <electronic/Everything.h>
#include <electronic/VanDerWaalsD3_data.h>


VanDerWaalsD3::VanDerWaalsD3(const Everything& e) : VanDerWaals(e)
{
	logPrintf("\nInitializing DFT-D3 calculator:\n");
	
	//Get parameters for exchange-correlation functional
	string xcName = e.exCorr.getName();
	D3::setXCscale(xcName, s6, sr6, s8, sr8);
	logPrintf("\tParameters set for %s functional\n", xcName.c_str());
	logPrintf("\ts6: %6.3lf  s_r6: %6.3lf\n", s6, sr6);
	logPrintf("\ts8: %6.3lf  s_r8: %6.3lf\n", s8, sr8);

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

