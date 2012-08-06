/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <core/DataIO.h>
#include <core/GridInfo.h>
#include <core/Operators.h>
#include <cblas.h>
#include <string.h>
#include <algorithm>


void saveDX(const DataRptr& X, const char* filenamePrefix)
{	char filename[256];
	sprintf(filename, "%s.bin", filenamePrefix);
	saveRawBinary(X, filename);

	sprintf(filename, "%s.dx", filenamePrefix);

	const GridInfo& g = X->gInfo;
	vector3<> Rcenter = 0.5*(g.R.column(0) + g.R.column(1) + g.R.column(2));
	vector3<> h0 = g.R.column(0)/g.S[0];
	vector3<> h1 = g.R.column(1)/g.S[1];
	vector3<> h2 = g.R.column(2)/g.S[2];

	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename)
	fprintf(fp, "object 1 class gridpositions counts %d %d %d\n", g.S[0], g.S[1], g.S[2]);
	fprintf(fp, "origin %lf %lf %lf\n", -Rcenter[0], -Rcenter[1], -Rcenter[2]);
	fprintf(fp, "delta %lf %lf %lf\n", h0[0], h0[1], h0[2]);
	fprintf(fp, "delta %lf %lf %lf\n", h1[0], h1[1], h1[2]);
	fprintf(fp, "delta %lf %lf %lf\n", h2[0], h2[1], h2[2]);
	fprintf(fp, "\n");

	fprintf(fp, "object 2 class gridconnections counts %d %d %d\n", g.S[0], g.S[1], g.S[2]);
	fprintf(fp, "attribute \"element type\" string \"cubes\"\n");
	fprintf(fp, "attribute \"ref\" string \"positions\"\n");
	fprintf(fp, "\n");

	fprintf(fp, "object 3 class array type double rank 0 items %d\n", g.nr);
	fprintf(fp, "ieee data file %s.bin\n", filenamePrefix);
	fprintf(fp, "\n");

	fprintf(fp, "object \"field\" class field\n");
	fprintf(fp, "component \"positions\" value 1\n");
	fprintf(fp, "component \"connections\" value 2\n");
	fprintf(fp, "component \"data\" value 3\n");
	fprintf(fp, "end\n");
	fprintf(fp, "\n");

	fclose(fp);
}


void saveSphericalized(const DataRptr* dataR, int nColumns, const char* filename, double drFac, vector3<>* centerPtr)
{	const GridInfo& g = dataR[0]->gInfo;

	//Determine the center for sphericalizaion:
	vector3<> center;
	if(centerPtr) center = *centerPtr;
	else center = 0.5*(g.R.column(0) + g.R.column(1) + g.R.column(2));

	//Determine the maximum radius from center to edges of the box:
	double rMax = 0.0;
	for(int i0=0; i0<2; i0++) for(int i1=0; i1<2; i1++) for(int i2=0; i2<2; i2++) //loop over vertices of parellelopiped
	{	double distance = (g.R.column(0)*i0 + g.R.column(1)*i1 + g.R.column(2)*i2 - center).length();
		if(distance > rMax) rMax = distance;
	}

	//Determine the diameter of the parellelopiped:
	double dr = 0.0;
	const vector3<>* h = g.h;
	vector3<> hCenter = 0.5*(h[0]+h[1]+h[2]);
	for(int i0=0; i0<2; i0++) for(int i1=0; i1<2; i1++) for(int i2=0; i2<2; i2++) //loop over vertices of parellelopiped
	{	double distance = (h[0]*i0 + h[1]*i1 + h[2]*i2 - hCenter).length();
		if(distance > dr) dr = distance;
	}
	dr *= (2.0*drFac);
	double drInv = 1.0/dr;

	//Sphericalize!
	int nRadial = int(ceil(rMax/dr));
	double** mean = new double*[nColumns];
	double* weight = new double[nRadial];
	for(int c=0; c<nColumns; c++)
	{	mean[c] = new double[nRadial];
		memset(weight, 0, nRadial*sizeof(double));
		memset(mean[c], 0, nRadial*sizeof(double));
		double* curData = dataR[c]->data();
		//Loop over all the lattice points and accumulate the above sums:
		int iStart=0, iStop=g.nr; const vector3<int> &S = g.S;
		THREAD_rLoop //not actually threaded, but convenient nonetheless
		(	vector3<> ri = iv[0]*h[0] + iv[1]*h[1] + iv[2]*h[2];
			double rRel = (ri - center).length() * drInv;
			int iRadial = int(floor(rRel));
			double wRight = (pow(rRel,2) - pow(iRadial,2))/(2*iRadial+1);
			double wLeft = 1.0 - wRight;
			if(wLeft && iRadial<nRadial)
			{	weight[iRadial] += wLeft;
				mean[c][iRadial] += wLeft * curData[i];
			}
			if(wRight && iRadial+1<nRadial)
			{	weight[iRadial+1] += wRight;
				mean[c][iRadial+1] += wRight * curData[i];
			}
		)
		//Calculate mean
		for(int i=0; i<nRadial; i++) mean[c][i] /= weight[i];
	}

	//Output data:
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename)
	for(int i=0; i<nRadial; i++)
	{	fprintf(fp, "%le", i*dr);
		for(int c=0; c<nColumns; c++) fprintf(fp, "\t%le", mean[c][i]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	//Cleanup:
	for(int c=0; c<nColumns; c++) delete[] mean[c];
	delete[] mean;
	delete[] weight;
}

void saveSphericalized(const DataGptr* dataG, int nColumns, const char* filename, double dGFac)
{	const GridInfo& g = dataG[0]->gInfo;
	int iStart=0, iStop=g.nG; const vector3<int>& S=g.S; const matrix3<>& GGT=g.GGT;

	double Gmax = 0.0;
	THREAD_halfGspaceLoop(
		double Gsq = GGT.metric_length_squared(iG);
		if(Gsq>Gmax) Gmax = Gsq;
	)
	Gmax = sqrt(Gmax);

	double dG = 0.0;
	for(iG[0]=-1; iG[0]<=1; iG[0]+=2) for(iG[1]=-1; iG[1]<=1; iG[1]+=2) for(iG[2]=-1; iG[2]<=1; iG[2]+=2) //loop over vertices of fft box
	{	double lengthG = sqrt(GGT.metric_length_squared(iG));
		if(lengthG > dG) dG = lengthG;
	}
	dG *= dGFac;
	double dGinv = 1.0/dG;

	//Sphericalize!
	int nRadial = int(ceil(Gmax/dG));
	double** mean = new double*[nColumns];
	double* weight = new double[nRadial];
	for(int c=0; c<nColumns; c++)
	{	mean[c] = new double[nRadial];
		memset(weight, 0, nRadial*sizeof(double));
		memset(mean[c], 0, nRadial*sizeof(double));
		complex* curData = dataG[c]->data();
		//Loop over all the reciprocal lattice points and accumulate the above sums:
		THREAD_halfGspaceLoop(
			double Grel = sqrt(GGT.metric_length_squared(iG)) * dGinv;
			int iRadial = int(floor(Grel));
			double wRight = (pow(Grel,2) - pow(iRadial,2))/(2*iRadial+1);
			double wLeft = 1.0 - wRight;
			if(wLeft && iRadial<nRadial)
			{	weight[iRadial] += wLeft;
				mean[c][iRadial] += wLeft * abs(curData[i]);
			}
			if(wRight && iRadial+1<nRadial)
			{	weight[iRadial+1] += wRight;
				mean[c][iRadial+1] += wRight * abs(curData[i]);
			}
		)
		//Calculate mean
		for(int i=0; i<nRadial; i++) mean[c][i] /= weight[i];
	}

	//Output data:
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename)
	for(int i=0; i<nRadial; i++)
	{	fprintf(fp, "%le", i*dG);
		for(int c=0; c<nColumns; c++) fprintf(fp, "\t%le", mean[c][i]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	//Cleanup:
	for(int c=0; c<nColumns; c++) delete[] mean[c];
	delete[] mean;
	delete[] weight;
}
