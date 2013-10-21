/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <core1D/Operators.h>
#include <core/BlasExtra.h>
#include <core/Random.h>

//------------------------------ Linear Unary operators ------------------------------

//Elementwise multiply on a ManagedMemory object (used by many operators below)
inline void dmul(const std::vector<double>& d, ManagedMemory& Y)
{	assert(Y);
	assert(d.size() == Y.nData());
	eblas_dmul(Y.nData(), d.data(),1, Y.data(),1);
}

//Elementwise divide on a ManagedMemory object (used by many operators below)
inline void ddiv(const std::vector<double>& d, ManagedMemory& Y)
{	assert(Y);
	assert(d.size() == Y.nData());
	eblas_ddiv(Y.nData(), d.data(),1, Y.data(),1);
}

//Row-major Dense-Matrix multiply (Y = A.X) on ManagedMemory objects (used by transform operators below)
inline void dgemv(bool transpose, const std::vector<double>& A, const ManagedMemory& X, ManagedMemory& Y)
{	assert(X);
	assert(Y);
	assert(A.size() == X.nData() * Y.nData());
	int M = transpose ? X.nData() : Y.nData();
	int N = transpose ? Y.nData() : X.nData();
	cblas_dgemv(CblasRowMajor, (transpose ? CblasTrans : CblasNoTrans),
		M, N, 1., A.data(), M,  X.data(),1, 0., Y.data(),1);
}

ScalarFieldTilde O(const ScalarFieldTilde& Y)
{	ScalarFieldTilde tmp(Y);
	return O((ScalarFieldTilde&&)tmp);
}
ScalarFieldTilde O(ScalarFieldTilde&& Y)
{	dmul(Y.gInfo->wTilde, Y);
	return Y;
}
ScalarFieldTilde Oinv(const ScalarFieldTilde& Y)
{	ScalarFieldTilde tmp(Y);
	return Oinv((ScalarFieldTilde&&)tmp);
}
ScalarFieldTilde Oinv(ScalarFieldTilde&& Y)
{	ddiv(Y.gInfo->wTilde, Y);
	return Y;
}

ScalarField I(const ScalarFieldTilde& Xtilde)
{	assert(Xtilde);
	const GridInfo& gInfo = *(Xtilde.gInfo);
	ScalarField X(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	ScalarFieldTilde tmp(Xtilde);
			dmul(gInfo.wTilde, tmp); //premultiply by basis weights
			dgemv(false, gInfo.matI, tmp, X); //multiply by matI
			break;
		}
		case GridInfo::Planar:
		{	fftw_execute_r2r(gInfo.planPlanarI, (double*)Xtilde.data(), X.data());
			X *= (1./gInfo.rMax);
			break;
		}
	}
	return X;
}

ScalarField ID(const ScalarFieldTilde& Xtilde)
{	assert(Xtilde);
	const GridInfo& gInfo = *(Xtilde.gInfo);
	ScalarField X(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	ScalarFieldTilde tmp(Xtilde);
			dmul(gInfo.wTilde, tmp); //premultiply by basis weights
			dgemv(false, gInfo.matID, tmp, X); //multiply by matID
			break;
		}
		case GridInfo::Planar:
		{	double normFac = (1./gInfo.rMax);
			//Gradient of cos Gr (obtained by doing the corresponding DST)
			//(need to shift down one point in G, and set final element to 0)
			ScalarFieldTilde tmp(&gInfo); double* tmpData = tmp.data();
			const double* XtildeData = Xtilde.data();
			for(int i=1; i<gInfo.S; i++)
				tmpData[i-1] = XtildeData[i] * (-normFac * gInfo.G[i]);
			tmpData[gInfo.S-1] = 0.;
			fftw_execute_r2r(gInfo.planPlanarID, (double*)tmp.data(), X.data());
			break;
		}
	}
	return X;
}

ScalarField IDD(const ScalarFieldTilde& Xtilde)
{	assert(Xtilde);
	const GridInfo& gInfo = *(Xtilde.gInfo);
	ScalarField X(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	ScalarFieldTilde tmp(Xtilde);
			dmul(gInfo.wTilde, tmp); //premultiply by basis weights
			dgemv(false, gInfo.matIDD, tmp, X); //multiply by matIDD
			break;
		}
		case GridInfo::Planar:
		{	double normFac = (1./gInfo.rMax);
			ScalarFieldTilde tmp(Xtilde); double* tmpData = tmp.data();
			for(int i=0; i<gInfo.S; i++)
				tmpData[i] *= (-normFac * gInfo.G[i]*gInfo.G[i]);
			fftw_execute_r2r(gInfo.planPlanarI, (double*)tmp.data(), X.data());
			break;
		}
	}
	return X;
}

ScalarField Jdag(const ScalarFieldTilde& Xtilde)
{	assert(Xtilde);
	const GridInfo& gInfo = *(Xtilde.gInfo);
	ScalarField X(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	dgemv(false, gInfo.matI, Xtilde, X); //multiply by matI
			dmul(gInfo.w, X); //postmultiply by quadrature weights
			break;
		}
		case GridInfo::Planar:
		{	ScalarFieldTilde tmp(Xtilde);
			tmp.data()[0] *= 2;
			fftw_execute_r2r(gInfo.planPlanarI, (double*)tmp.data(), X.data());
			X *= (0.5 * gInfo.rMax / gInfo.S);
			break;
		}
	}
	return X;
}



ScalarFieldTilde J(const ScalarField& X)
{	assert(X);
	const GridInfo& gInfo = *(X.gInfo);
	ScalarFieldTilde Xtilde(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	ScalarField tmp(X);
			dmul(gInfo.w, tmp); //premultiply by quadrature weights
			dgemv(true, gInfo.matI, tmp, Xtilde); //multiply by transpose(matI)
			break;
		}
		case GridInfo::Planar:
		{	fftw_execute_r2r(gInfo.planPlanarIdag, (double*)X.data(), Xtilde.data());
			Xtilde *= (0.5 * gInfo.rMax / gInfo.S);
			break;
		}
	}
	return Xtilde;
}

ScalarFieldTilde Idag(const ScalarField& X)
{	assert(X);
	const GridInfo& gInfo = *(X.gInfo);
	ScalarFieldTilde Xtilde(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	dgemv(true, gInfo.matI, X, Xtilde); //multiply by transpose(matI)
			dmul(gInfo.wTilde, Xtilde); //postmultiply by basis weights
			break;
		}
		case GridInfo::Planar:
		{	fftw_execute_r2r(gInfo.planPlanarIdag, (double*)X.data(), Xtilde.data());
			Xtilde *= (1./gInfo.rMax);
			Xtilde.data()[0] *= 0.5;
			break;
		}
	}
	return Xtilde;
}

ScalarFieldTilde IDdag(const ScalarField& X)
{	assert(X);
	const GridInfo& gInfo = *(X.gInfo);
	ScalarFieldTilde Xtilde(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	dgemv(true, gInfo.matID, X, Xtilde); //multiply by transpose(matID)
			dmul(gInfo.wTilde, Xtilde); //postmultiply by basis weights
			break;
		}
		case GridInfo::Planar:
		{	double normFac = (1./gInfo.rMax);
			//(transpose of) gradient of cos kr (appropriate DST)
			double* XtildeData = Xtilde.data();
			fftw_execute_r2r(gInfo.planPlanarIDdag, (double*)X.data(), XtildeData);
			//(need to shift up one point in k in the result)
			for(int i=gInfo.S-1; i>0; i--)
				XtildeData[i] = XtildeData[i-1] * (-normFac * gInfo.G[i]);
			XtildeData[0] = 0.;
			break;
		}
	}
	return Xtilde;
}

ScalarFieldTilde IDDdag(const ScalarField& X)
{	assert(X);
	const GridInfo& gInfo = *(X.gInfo);
	ScalarFieldTilde Xtilde(&gInfo);
	switch(gInfo.coord)
	{
		case GridInfo::Spherical:
		case GridInfo::Cylindrical:
		{	dgemv(true, gInfo.matIDD, X, Xtilde); //multiply by transpose(matIDD)
			dmul(gInfo.wTilde, Xtilde); //postmultiply by basis weights
			break;
		}
		case GridInfo::Planar:
		{	double normFac = (1./gInfo.rMax);
			double* XtildeData = Xtilde.data();
			fftw_execute_r2r(gInfo.planPlanarIdag, (double*)X.data(), XtildeData);
			for(int i=0; i<gInfo.S; i++)
				XtildeData[i] *= (-normFac * gInfo.G[i]*gInfo.G[i]);
			break;
		}
	}
	return Xtilde;
}



ScalarField JdagOJ(const ScalarField& Y)
{	ScalarField tmp(Y);
	return JdagOJ((ScalarField&&)tmp);
}
ScalarField JdagOJ(ScalarField&& Y)
{	dmul(Y.gInfo->w, Y);
	return Y;
}
ScalarField DiagJdagOJ1(const ScalarField& Y)
{	ScalarField tmp(Y);
	return DiagJdagOJ1((ScalarField&&)tmp);
}
ScalarField DiagJdagOJ1(ScalarField&& Y)
{	dmul(Y.gInfo->w, Y);
	return Y;
}
ScalarField DiagJdagOJ1(double s, const GridInfo& gInfo)
{	ScalarField ret(&gInfo);
	ret.zero();
	eblas_daxpy(gInfo.S, s, gInfo.w.data(),1, ret.data(),1);
	return ret;
}
ScalarField DiagJdagOJ1inv(const ScalarField& Y)
{	ScalarField tmp(Y);
	return DiagJdagOJ1inv((ScalarField&&)tmp);
}
ScalarField DiagJdagOJ1inv(ScalarField&& Y)
{	ddiv(Y.gInfo->w, Y);
	return Y;
}

ScalarFieldTilde L(const ScalarFieldTilde& Y)
{	ScalarFieldTilde tmp(Y);
	return L((ScalarFieldTilde&&)tmp);
}
ScalarFieldTilde L(ScalarFieldTilde&& Y)
{	assert(Y);
	const GridInfo& gInfo = *(Y.gInfo);
	double* Ydata = Y.data();
	for(int i=0; i<gInfo.S; i++)
		Ydata[i] *= (-gInfo.wTilde[i] * gInfo.G[i] * gInfo.G[i]);
	return Y;
}
ScalarFieldTilde Linv(const ScalarFieldTilde& Y)
{	ScalarFieldTilde tmp(Y);
	return Linv((ScalarFieldTilde&&)tmp);
}
ScalarFieldTilde Linv(ScalarFieldTilde&& Y)
{	assert(Y);
	const GridInfo& gInfo = *(Y.gInfo);
	double* Ydata = Y.data();
	Ydata[0] = 0.;
	for(int i=1; i<gInfo.S; i++)
		Ydata[i] /= (-gInfo.wTilde[i] * gInfo.G[i] * gInfo.G[i]);
	return Y;
}

//------------------------------ Nonlinear Unary operators ------------------------------

ScalarField exp(const ScalarField& Y)
{	assert(Y);
	ScalarField ret(Y.gInfo);
	double* retData = ret.data();
	const double* Ydata = Y.data();
	for(unsigned i=0; i<Y.nData(); i++)
		retData[i] = exp(Ydata[i]);
	return ret;
}
ScalarField log(const ScalarField& Y)
{	assert(Y);
	ScalarField ret(Y.gInfo);
	double* retData = ret.data();
	const double* Ydata = Y.data();
	for(unsigned i=0; i<Y.nData(); i++)
		retData[i] = log(Ydata[i]);
	return ret;
}
ScalarField sqrt(const ScalarField& Y)
{	assert(Y);
	ScalarField ret(Y.gInfo);
	double* retData = ret.data();
	const double* Ydata = Y.data();
	for(unsigned i=0; i<Y.nData(); i++)
		retData[i] = sqrt(Ydata[i]);
	return ret;
}
ScalarField inv(const ScalarField& Y)
{	assert(Y);
	ScalarField ret(Y.gInfo);
	double* retData = ret.data();
	const double* Ydata = Y.data();
	for(unsigned i=0; i<Y.nData(); i++)
		retData[i] = 1./Ydata[i];
	return ret;
}
ScalarField pow(const ScalarField& Y, double alpha)
{	assert(Y);
	ScalarField ret(Y.gInfo);
	double* retData = ret.data();
	const double* Ydata = Y.data();
	for(unsigned i=0; i<Y.nData(); i++)
		retData[i] = pow(Ydata[i], alpha);
	return ret;
}


//------------------------------ Multiplication operators ------------------------------

//scale and unary-
ScalarField& operator*=(ScalarField& X, double s) { scale(s, X); return X; }
scaled<ScalarField> operator*(double s, const ScalarField &Y) { return scaled<ScalarField>(Y, s); }
scaled<ScalarField> operator*(const ScalarField &Y, double s) { return scaled<ScalarField>(Y, s); }
scaled<ScalarField> operator-(const ScalarField &Y) { return scaled<ScalarField>(Y, -1); }
ScalarFieldTilde& operator*=(ScalarFieldTilde& X, double s) { scale(s, X); return X; }
scaled<ScalarFieldTilde> operator*(double s, const ScalarFieldTilde &Y) { return scaled<ScalarFieldTilde>(Y, s); }
scaled<ScalarFieldTilde> operator*(const ScalarFieldTilde &Y, double s) { return scaled<ScalarFieldTilde>(Y, s); }
scaled<ScalarFieldTilde> operator-(const ScalarFieldTilde &Y) { return scaled<ScalarFieldTilde>(Y, -1); }

//Convolution:
ScalarFieldTilde operator*(const SphericalKernel& K, const ScalarFieldTilde& Y)
{	ScalarFieldTilde tmp(Y);
	return K * ((ScalarFieldTilde&&)tmp);
}
ScalarFieldTilde operator*(const SphericalKernel& K, ScalarFieldTilde&& Y)
{	dmul(K, Y);
	return Y;
}

//Elementwise multiply
ScalarField DiagScalarField::operator*(const ScalarField& X) const
{	ScalarField Xcopy(X);
	return (*this) * ((ScalarField&&)Xcopy);
}
ScalarField DiagScalarField::operator*(ScalarField&& X) const
{	assert(X.gInfo == d.gInfo);
	eblas_dmul(X.nData(), d.data(),1, X.data(),1);
	return X;
}
DiagScalarField Diag(const ScalarField& X)
{	assert(X);
	return DiagScalarField(X);
}

//------------------------------ Linear combine operators ------------------------------

ScalarField& operator+=(ScalarField& Y, const ScalarField &X) { if(Y) { if(X) axpy(1.0, X, Y); } else Y=X; return Y; }
ScalarField& operator-=(ScalarField& Y, const ScalarField &X) { if(Y) { if(X) axpy(-1.0, X, Y); } else Y=-X; return Y; }
ScalarField operator+(const ScalarField &Y1, const ScalarField &Y2) { ScalarField Ysum(Y1); Ysum += Y2; return Ysum; }
ScalarField operator-(const ScalarField &Y1,const ScalarField &Y2) { ScalarField Ydiff(Y1); Ydiff -= Y2; return Ydiff; }
ScalarFieldTilde& operator+=(ScalarFieldTilde& Y, const ScalarFieldTilde &X) { if(Y) { if(X) axpy(1.0, X, Y); } else Y=X; return Y; }
ScalarFieldTilde& operator-=(ScalarFieldTilde& Y, const ScalarFieldTilde &X) { if(Y) { if(X) axpy(-1.0, X, Y); } else Y=-X; return Y; }
ScalarFieldTilde operator+(const ScalarFieldTilde &Y1, const ScalarFieldTilde &Y2) { ScalarFieldTilde Ysum(Y1); Ysum += Y2; return Ysum; }
ScalarFieldTilde operator-(const ScalarFieldTilde &Y1,const ScalarFieldTilde &Y2) { ScalarFieldTilde Ydiff(Y1); Ydiff -= Y2; return Ydiff; }

//Real-space scalar additions:
ScalarField& operator+=(ScalarField& Y, double a)
{	eblas_daxpy(Y.nData(), 1., &a, 0, Y.data(), 1);
	return Y;
}
ScalarField operator+(double a, const ScalarField& Y) { ScalarField X(Y); return X += a; }
ScalarField operator+(const ScalarField& Y, double a) { ScalarField X(Y); return X += a; }
ScalarField& operator-=(ScalarField& Y, double a) { return Y += -a; }
ScalarField operator-(double a, const ScalarField& Y) { ScalarField X(Y); return (X -= a) *= (-1.); }
ScalarField operator-(const ScalarField& Y, double a) { ScalarField X(Y); return X -= a; }

double integral(const ScalarField& x)
{	assert(x);
	const GridInfo& gInfo = *(x.gInfo);
	return eblas_ddot(gInfo.S, gInfo.w.data(),1, x.data(),1);
}

double integral(const ScalarFieldTilde& xTilde)
{	assert(xTilde);
	return xTilde.data()[0]; //Just the G=0 component for all plane-wave related bases
}


//------------------------- Initialization utilties ---------------------------------

void nullToZero(ScalarField& X, const GridInfo& gInfo)
{	if(!X) X.init(&gInfo);
}

void initRandom(ScalarField& X, double cap)
{	double* Xdata = X.data();
	for(unsigned i=0; i<X.nData(); i++)
		Xdata[i] = Random::normal(0, 1, cap);
}

void initRandomFlat(ScalarField& X)
{	double* Xdata = X.data();
	for(unsigned i=0; i<X.nData(); i++)
		Xdata[i] = Random::uniform();
}
