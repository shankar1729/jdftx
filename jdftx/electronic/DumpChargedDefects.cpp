#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/Dump_internal.h>
#include <core/WignerSeitz.h>
#include <core/Operators.h>
#include <core/DataIO.h>

//-------------------------- Slab epsilon ----------------------------------

inline void planarAvg_sub(size_t iStart, size_t iStop, const vector3<int>& S, int iDir, complex* data)
{	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	THREAD_halfGspaceLoop( if(iG[jDir] || iG[kDir]) data[i] = 0.; )
}

inline void fixBoundary_sub(size_t iStart, size_t iStop, const vector3<int>& S, int iDir, int iBoundary, double* eps)
{	iBoundary = positiveRemainder(iBoundary, S[iDir]);
	THREAD_rLoop(
		int dist = iv[iDir] - iBoundary;
		if(dist*2 > S[iDir]) dist -= S[iDir];
		if(dist*2 < -S[iDir]) dist += S[iDir];
		if(abs(dist)<=2) eps[i] = 1.;
	)
}

void SlabEpsilon::dump(const Everything& e, DataRptr d_tot) const
{	string fname = e.dump.getFilename("slabEpsilon");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	//Read reference Dtot:
	DataRptr d_totRef(DataR::alloc(e.gInfo));
	loadRawBinary(d_totRef, dtotFname.c_str());
	//Calculate inverse of epsilon:
	int iDir = e.coulombParams.iDir;
	vector3<> zHat = e.gInfo.R.column(iDir);
	zHat *= 1./zHat.length(); //unit vector along slab normal
	double dE = dot(zHat, e.coulombParams.Efield - Efield);
	if(!dE) die("\nThe applied electric fields in the reference and present calculations are equal.\n");
	//--- calculate field using central-difference derivative:
	DataGptr tPlus(DataG::alloc(e.gInfo)), tMinus(DataG::alloc(e.gInfo));
	initTranslation(tPlus, e.gInfo.h[iDir]);
	initTranslation(tMinus, -e.gInfo.h[iDir]);
	double h = e.gInfo.h[iDir].length();
	DataGptr epsInvTilde = (1./(dE * 2.*h)) * (tPlus - tMinus) * J(d_tot - d_totRef);
	threadLaunch(planarAvg_sub, e.gInfo.nG, e.gInfo.S, iDir, epsInvTilde->data()); //planar average
	//Fix values of epsilon near truncation boundary to 1:
	DataRptr epsInv = I(epsInvTilde);
	threadLaunch(fixBoundary_sub, e.gInfo.nr, e.gInfo.S, iDir, e.coulomb->ivCenter[iDir] + e.gInfo.S[iDir]/2, epsInv->data());
	//Apply smoothing:
	epsInv = I(gaussConvolve(J(epsInv), sigma));
	//Write file:
	FILE* fp = fopen(fname.c_str(), "w");
	fprintf(fp, "#distance[bohr]  epsilon\n");
	vector3<int> iR;
	double* epsInvData = epsInv->data();
	for(iR[iDir]=0; iR[iDir]<e.gInfo.S[iDir]; iR[iDir]++)
		fprintf(fp, "%lf %lf\n", h*iR[iDir], 1./epsInvData[e.gInfo.fullRindex(iR)]);
	fclose(fp);
	logPrintf("done\n"); logFlush();
}

//-------------------------- Charged defects ----------------------------------


void ChargedDefect::dump(const Everything& e, DataRptr d_tot) const
{	logPrintf("Calculating charged defect correction:\n"); logFlush();
	
	//Construct model charge on plane-wave grid:
	DataGptr rhoModel(DataG::alloc(e.gInfo));
	initTranslation(rhoModel, e.gInfo.R * pos);
	rhoModel = gaussConvolve((q/e.gInfo.detR)*rhoModel, sigma);
	
	//Calculate isolated and periodic self-energy (and potential) of model charge
	DataGptr dModel; double Emodel=0., EmodelIsolated=0.;
	switch(e.coulombParams.geometry)
	{	case CoulombParams::Periodic: //Bulk defect
		{	//Periodic potential and energy:
			dModel = (*e.coulomb)(rhoModel) * (1./bulkEps); //assuming uniform dielectric
			Emodel = 0.5*dot(rhoModel, O(dModel));
			//Isolated energy:
			EmodelIsolated = 0.5*q*q/(sigma*sqrt(M_PI)*bulkEps); //self energy of Gaussian
			break;
		}
		case CoulombParams::Slab: //Surface defect
		{	die("Not yet implemented.\n");
			break;
		}
		default: die("Coulomb-interaction geometry must be either slab or periodic for charged-defect correction.\n");
	}
	logPrintf("\tEmodelIsolated: %.8lf\n", EmodelIsolated);
	logPrintf("\tEmodelPeriodic: %.8lf\n", Emodel);
	
	//Calculate alignment potential:
	//--- read reference Dtot and calculate electrostatic potentials:
	DataRptr d_totRef(DataR::alloc(e.gInfo));
	loadRawBinary(d_totRef, dtotFname.c_str());
	DataRptr Vdft = d_tot - d_totRef; //electrostatic potential of defect from DFT
	DataRptr Vmodel = I(dModel);  //electrostatic potential from model
	//--- radially histogram:
	DataRptr Varr[2] = { Vdft, Vmodel };
	vector3<> rCenter = e.gInfo.R*pos;
	std::vector< std::vector<double> > hist = sphericalize(Varr, 2, 1., &rCenter);
	const std::vector<double>& rRadial = hist[0];
	const std::vector<double>& VdftRadial = hist[1];
	const std::vector<double>& VmodelRadial = hist[2];
	const std::vector<double>& wRadial = hist[3];
	//--- calculate alignment scalar:
	double deltaV = 0.; double wSum = 0.;
	for(size_t i=0; i<rRadial.size(); i++)
	{	double r = rRadial[i];
		double w = wRadial[i] * 0.5*erfc((rMin-r)/rSigma);
		deltaV += (VmodelRadial[i] - VdftRadial[i])*w;
		wSum += w;
	}
	deltaV /= wSum;
	logPrintf("\tDeltaV:         %.8lf (average over %lg grid points)\n", deltaV, wSum);
	logPrintf("\tNet correction: %.8lf (= EmodelIsolated - EmodelPeriodic + q deltaV)\n", EmodelIsolated - Emodel + q*deltaV);
	//--- write alignment data
	string fname = e.dump.getFilename("chargedDefectDeltaV");
	logPrintf("\tWriting %s (plot to check DeltaV manually) ... ", fname.c_str()); logFlush();
	FILE* fp = fopen(fname.c_str(), "w");
	if(!fp) die("Error opening %s for writing.\n", fname.c_str())
	fprintf(fp, "#r DeltaV Vmodel Vdft\n");
	for(size_t i=0; i<rRadial.size(); i++)
		fprintf(fp, "%lf %le %le %le\n", rRadial[i], VmodelRadial[i]-VdftRadial[i], VmodelRadial[i], VdftRadial[i]);
	fclose(fp);
	logPrintf("done.\n");
}
