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

#include <core/Random.h>
#include <random>

namespace Random
{
	std::mt19937_64 generator;
	std::normal_distribution<double> normdist;
	std::uniform_real_distribution<double> uniformDist;
	
	void seed(int i)
	{	generator.seed(i);
	}
	
	double uniform(double start, double end)
	{	return start + (end-start)*uniformDist(generator);
	}

	double normal(double mean, double sigma, double cap)
	{	double r = sigma*normdist(generator);
		if(cap>0.0) while(fabs(r)>cap) r = sigma*normdist(generator); //regenerate till cap is satisfied
		return mean + r;
	}

	complex normalComplex(double sigma)
	{	return sigma * complex(normdist(generator), normdist(generator));
	}
}
