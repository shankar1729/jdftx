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

#include <fluid/SO3quad.h>
#include <fluid/Euler.h>
#include <string.h>

SO3quad::SO3quad(const S2quad& s2quad, int Zn)
{	setup(s2quad, Zn);
}

SO3quad::SO3quad(S2quadType type, int Zn, unsigned nBeta, unsigned nAlpha, unsigned nGamma)
{	S2quad* s2quad = 0;
	switch(type)
	{	case QuadEuler:        s2quad = new EulerProduct(nBeta, nAlpha, nGamma); break;
		case QuadTetrahedron:  s2quad = new Tetrahedron(); break;
		case QuadOctahedron:   s2quad = new Octahedron(); break;
		case QuadIcosahedron:  s2quad = new Icosahedron(); break;
		case Quad7design_24:   s2quad = new S2_7design_24(); break;
		case Quad8design_36:   s2quad = new S2_8design_36(); break;
		case Quad9design_48:   s2quad = new S2_9design_48(); break;
		case Quad10design_60:  s2quad = new S2_10design_60(); break;
		case Quad11design_70:  s2quad = new S2_11design_70(); break;
		case Quad12design_84:  s2quad = new S2_12design_84(); break;
		case Quad13design_94:  s2quad = new S2_13design_94(); break;
		case Quad14design_108: s2quad = new S2_14design_108(); break;
		case Quad15design_120: s2quad = new S2_15design_120(); break;
		case Quad16design_144: s2quad = new S2_16design_144(); break;
		case Quad17design_156: s2quad = new S2_17design_156(); break;
		case Quad18design_180: s2quad = new S2_18design_180(); break;
		case Quad19design_204: s2quad = new S2_19design_204(); break;
		case Quad20design_216: s2quad = new S2_20design_216(); break;
		case Quad21design_240: s2quad = new S2_21design_240(); break;
	}
	setup(*s2quad, Zn);
	delete s2quad;
}


void SO3quad::print()
{	for(unsigned o=0; o<eulerS2.size(); o++)
		logPrintf("Orientations %3u - %3u with weight %lf:  [ %+.3lf %+.3lf %+.3lf+2n/%d ] x pi\n",
			nS1byZn*o+1, nS1byZn*(o+1), weightS2[o],
			eulerS2[o][0]/M_PI, eulerS2[o][1]/M_PI, eulerS2[o][2]/M_PI, nS1);
}

int SO3quad::nOrientations() const
{	return weightS2.size() * nS1byZn;
}

vector3<> SO3quad::euler(int iOrientation) const
{	vector3<> ret = eulerS2[iOrientation/nS1byZn]; //get the base euler angle
	ret[2] += (2*M_PI/nS1) * (iOrientation%nS1byZn); //add the S1 phase
	return ret;
}

double SO3quad::weight(int iOrientation) const
{	return weightS2[iOrientation/nS1byZn];
}

void SO3quad::setup(const S2quad& s2quad, int Zn)
{
	nS1byZn = ceil(s2quad.nS1()*1./Zn);
	nS1 = nS1byZn * Zn; //round up to nearest multiple of Zn
	eulerS2 = s2quad.euler; //copy over S2 quadrature nodes
	weightS2 = s2quad.weight; //copy over weights
	assert(eulerS2.size()==weightS2.size());
	
	logPrintf("Generating SO(3)/Z%d quadrature with %d nodes:\n", Zn, nOrientations());
	
	//Make sure the weights add up to 1 (including S1 multiplicity):
	int nWeights = eulerS2.size();
	double weightSum = 0.0;
	for(int i=0; i<nWeights; i++) weightSum += weightS2[i];
	for(int i=0; i<nWeights; i++) weightS2[i] *= 1./(weightSum * nS1byZn);
	logPrintf("\tFrom %d-point quadrature '%s' on S2 x %d uniform samples on S1/Z%d.\n",
		nWeights, s2quad.name().c_str(), nS1byZn, Zn);
	
	const int jMax = s2quad.jMax();
	if(!jMax) return;
	
	logPrintf("\tVerifying exactness to jMax = %d ... ", jMax); logFlush();
	double rmsErr=0., maxErr=0.; int nEquations=0;
	const double sumThreshold = 1e-6;
	for(int j=1; j<=jMax; j++)
	{	//Note only m2 = multiples of nS1 yield non-trivial equations for uniform S1 sampling
		//Consequently non-zero m2 gets checked only when nS1 <= jMax
		//(Amongst the currently implemented quadratures, this happens only for the Icosahedron group)
		int m2max = (j/nS1)*nS1;
		for(int m1=-j; m1<=j; m1++) for(int m2=-m2max; m2<=m2max; m2+=nS1)
		{	//Equations labelled by j, m1, m2: now loop over S2 weights
			double sum = 0.;
			for(int iWeight=0; iWeight<nWeights; iWeight++)
			{	double dj_m1m2 = wigner_d(j,m1,m2, eulerS2[iWeight][1]);
				double theta = m1 * eulerS2[iWeight][0] + m2 * eulerS2[iWeight][2];
				sum += nS1byZn * weightS2[iWeight] * dj_m1m2 * (m1>=0 ? cos(theta) : sin(theta));
			}
			if(fabs(sum) > sumThreshold)
				die("\n\tInexact at j=%d, m1=%d, m2=%d: |error|~%.0le > threshold~%.0le",
					j, m1, m2, fabs(sum), sumThreshold)
			rmsErr += sum*sum;
			if(fabs(sum)>maxErr) maxErr = fabs(sum);
			nEquations++;
		}
	}
	rmsErr = sqrt(rmsErr/nEquations);
	logPrintf("Done (MaxError~%.0le, RMSerror~%.0le.)\n", maxErr, rmsErr); logFlush();
}

