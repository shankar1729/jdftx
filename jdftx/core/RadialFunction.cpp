/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <core/RadialFunction.h>
#include <core/SphericalHarmonics.h>
#include <core/GpuUtil.h>
#include <core/Thread.h>

RadialFunctionG::RadialFunctionG() : dGinv(0), nCoeff(0),
#ifdef GPU_ENABLED
coeffGpu(0),
#endif
rFunc(0)
{
}

void RadialFunctionG::init(int l, int nSamples, double dG, const char* filename, double scale)
{	std::vector<double> samples(nSamples);
	FILE* fp = fopen(filename, "r");
	if(!fp) die("Could not open radial function file '%s' for reading.\n", filename);
	for(int i=0; i<nSamples; i++)
	{	if(fscanf(fp, "%*f %lf", &samples[i])<1)
			die("Error reading sample# %d of %d from radial function file '%s'.\n", i, nSamples, filename);
		samples[i] *= scale; //apply optional scale factor
	}
	fclose(fp);
	init(l, samples, dG);
}

void RadialFunctionG::init(int l, const std::vector<double>& samples, double dG)
{	set(QuinticSpline::getCoeff(samples, l%2==1), 1./dG);
}

void RadialFunctionG::set(const std::vector<double>& coeff, double dGinv)
{	this->coeff = coeff;
	this->dGinv = dGinv;
	nCoeff = coeff.size();
	#ifdef GPU_ENABLED
	cudaMalloc(&coeffGpu, sizeof(double)*nCoeff);
	cudaMemcpy(coeffGpu, coeff.data(), sizeof(double)*nCoeff, cudaMemcpyHostToDevice);
	gpuErrorCheck();
	#endif
}

void RadialFunctionG::updateGmax(int l, int nSamples)
{	if(!rFunc) return; //don't have the information necessary to update
	if(nSamples+4 <= nCoeff) return; //already have enough samples
	rFunc->transform(l, 1./dGinv, nSamples, *this);
}

void RadialFunctionG::free(bool rFuncDelete)
{	if(rFunc && rFuncDelete) delete rFunc;
	#ifdef GPU_ENABLED
	if(nCoeff) cudaFree(coeffGpu);
	#endif
}

RadialFunctionR::RadialFunctionR(int nSamples) : r(nSamples), dr(nSamples), f(nSamples)
{
}
RadialFunctionR::RadialFunctionR(const std::vector<double>& r, double dlogr)
: r(r), dr(r.size()), f(r.size())
{	for(unsigned i=0; i<r.size(); i++)
	{	dr[i] = r[i] * dlogr;
		if(i && fabs(log(r[i]/r[i-1]) - dlogr) > 1e-6)
			die("Logarithmic grid logPrintf(r[%d]/r[%d]) != %lg (specified dlogr)\n", i, i-1, dlogr);
	}
}
void RadialFunctionR::set(std::vector<double> r, std::vector<double> dr)
{	assert(r.size() >= this->f.size());
	assert(dr.size() >= this->f.size());
	this->r.assign(r.begin(), r.begin()+this->f.size());
	this->dr.assign(dr.begin(), dr.begin()+this->f.size());
}

void RadialFunctionR_initWeights_sub(int jStart, int jStop, int nSamples, const double* r, double* dr)
{	for(int j=jStart; j<jStop; j++)
	{	//Calculate cubic spline with natural B.C. that is 1 at j'th point and 0 at all others:
		std::vector<double> a(nSamples, 0.); //sub-diagonal (first entry ignored)
		std::vector<double> b(nSamples, 0.); //diagonal
		std::vector<double> c(nSamples, 0.); //super-diagonal (last entry ignored)
		std::vector<double> d(nSamples, 0.); //r.h.s
		std::vector<double> y(nSamples, 0.); y[j]=1; //values of spline y(r) at nodes
		for(int i=0; i<nSamples-1; i++)
		{	//Contribution from (r[i],r[i+1]) interval:
			double hInv = 1./(r[i+1]-r[i]);
			double rhs = 3*hInv*hInv*(y[i+1]-y[i]);
			b[i] += 2*hInv;   c[i] += hInv;     d[i] += rhs; //to i'th equation
			a[i+1] += hInv; b[i+1] += 2*hInv; d[i+1] += rhs; //to i+1'th equation
		}
		//Solve tridiagonal matrix equation to fully determine spline y(r):
		for(int i=1; i<nSamples; i++)
		{	double frac = a[i]/b[i-1];
			b[i] -= frac*c[i-1];
			d[i] -= frac*d[i-1];
		}
		std::vector<double> yp(nSamples); //values of dy/dr at each node
		yp[nSamples-1] = d[nSamples-1]/b[nSamples-1];
		for(int i=nSamples-2; i>=0; i--)
			yp[i] = (d[i] - c[i]*yp[i+1])/b[i];
		
		//Integrate spline over entire domain to get quadrature weight for current node:
		double r0 = r[0], y0 = y[0], yp0 = yp[0];
		double wj = (1./3)*(r0*r0*r0)*(y0-0.25*r0*yp0); //contribution from r=0 till first node
		for(int i=0; i<nSamples-1; i++)
		{	double r1 = r[i+1], y1 = y[i+1], yp1 = yp[i+1];
			//Contribution from (r[i],r[i+1]) interval:
			wj += (1./60) * (r1-r0) * ( r0*r1* 10*(y0+y1)
					+ r0*r0 * (16*y0 + 4*y1 + r0*(yp1-2*yp0) + r1*yp1)
					+ r1*r1 * (16*y1 + 4*y0 + r1*(yp0-2*yp1) + r0*yp0) );
			r0 = r1; y0 = y1; yp0 = yp1;
		}
		dr[j] = wj / (r[j]*r[j]);
	}
}

void RadialFunctionR::initWeights()
{	static StopWatch watch("initWeights"); watch.start();
	//Ensure r is valid (positive and in ascending order):
	int nSamples = r.size();
	assert(nSamples>=2);
	assert(r[0]>=0.);
	for(int i=0; i<nSamples-1; i++)
		assert(r[i+1]>r[i]);
	dr.resize(nSamples);
	if(!r[0]) r[0]=1e-6*r[1]; //avoid division by 0 (will make no difference at double prec)
	threadLaunch(RadialFunctionR_initWeights_sub, nSamples, nSamples, r.data(), dr.data());
	watch.stop();
}


// Evaluate spherical bessel transform of order l at wavevector G
double RadialFunctionR::transform(int l, double G) const
{	//Simpson's 1/3 rule on the log-grid:
	double sum = 0.0;
	int Nhlf = (f.size()-1)/2; //half the sample points
	for(int i=0; i<=2*Nhlf; i++)
	{	sum += f[i] * bessel_jl(l, G*r[i]) //sample value for spherical bessel transform
			* (dr[i] * r[i]*r[i]) //weight for radial spherical integration
			* ( (i==0 || i==2*Nhlf) ? 1 : 2*((i%2)+1) ); //simpson integration factor
	}
	return (4*M_PI/3) * sum;
}

static void RadialFunction_transform_sub(int iMineStart, int iMineStop, int iGstart, int l, double dG, const RadialFunctionR* rFunc, double* fTilde)
{	for(int iG=iGstart+iMineStart; iG<iGstart+iMineStop; iG++)
		fTilde[iG] = rFunc->transform(l, iG*dG);
}

// Initialize a uniform G radial function from the log-grid function
void RadialFunctionR::transform(int l, double dG, int nGrid, RadialFunctionG& func) const
{	static StopWatch watch("RadialFunctionR::transform"); watch.start();
	std::vector<double> fTilde(nGrid, 0.);
	int iGstart, iGstop; TaskDivision(nGrid, mpiWorld).myRange(iGstart, iGstop);
	int nGridMine = iGstop-iGstart;
	if(nGridMine)
		threadLaunch(RadialFunction_transform_sub, nGridMine, iGstart, l, dG, this, fTilde.data());
	mpiWorld->allReduceData(fTilde, MPIUtil::ReduceSum);
	func.free(this!=func.rFunc);
	func.init(l, fTilde, dG);
	if(this!=func.rFunc) func.rFunc = new RadialFunctionR(*this);
	watch.stop();
}

