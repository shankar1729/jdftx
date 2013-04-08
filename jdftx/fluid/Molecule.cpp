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

#include <fluid/ErfFMTweight.h>
#include <fluid/Molecule.h>
#include <core/Operators.h>


int get_nSites(std::vector<std::vector<vector3<>>>& PositionList)
{
	int counter = 0;
	for (unsigned i=0; i<PositionList.size(); i++){	counter += PositionList[i].size();	}
	return counter;
}

std::vector<Site> get_Site(std::vector<SiteProperties*>& PropList,std::vector<std::vector<vector3<>>>& PositionList)
{
	std::vector<Site> site;
	for (unsigned i=0; i<PropList.size(); i++)
	{
		for (unsigned j=0; j<PositionList[i].size(); j++) { site.push_back(Site(i, PropList[i], PositionList[i][j])); }
	}
	return site;	
}

Molecule::Molecule(std::vector<SiteProperties*>& PropList, std::vector<std::vector<vector3<>>>& PositionList, string name)
: name(name), site(get_Site(PropList,PositionList)), nSites(get_nSites(PositionList)), nIndices(PropList.size())
{
}

SiteProperties::SiteProperties(const GridInfo& gInfo, double sphereRadius, double sphereSigma, double chargeZ,
	RealKernel* chargeKernel, bool indepSite)
: sphereRadius(sphereRadius), sphereSigma(sphereSigma), chargeZ(chargeZ), chargeKernel(chargeKernel), indepSite(indepSite),
couplingZnuc(0), atomicNumber(0), couplingElecKernel(0)
{	
	if(sphereRadius)
	{	w0 = new RealKernel(gInfo);
		w1 = new RealKernel(gInfo);
		w2 = new RealKernel(gInfo);
		w3 = new RealKernel(gInfo);
		w1v = new RealKernel(gInfo);
		w2m = new RealKernel(gInfo);
		ErfFMTweight erfFMTweight(sphereRadius, sphereSigma);
		applyFuncGsq(gInfo, erfFMTweight, w0->data, w1->data, w2->data, w3->data, w1v->data, w2m->data);
		w0->set();
		w1->set();
		w2->set();
		w3->set();
		w1v->set();
		w2m->set();
	}
}

//probably will go when H2Osites gets replaced by SiteProperties
SiteProperties::SiteProperties(const GridInfo& gInfo, double sphereRadius, double sphereSigma,
	H2OSite& water_site, RealKernel* chargeKernel, bool indepSite)
: sphereRadius(sphereRadius), sphereSigma(sphereSigma), chargeZ(water_site.Z), chargeKernel(chargeKernel), indepSite(indepSite),
couplingZnuc(water_site.Znuc), atomicNumber(water_site.atomicNumber), couplingElecKernel(0), siteName(water_site.name), convCouplingWidth(water_site.CouplingWidth),
convCouplingSiteCharge(water_site.Z), convCouplingSiteModel(water_site.ccSiteModel) 
{
	if(sphereRadius)
	{	w0 = new RealKernel(gInfo);
		w1 = new RealKernel(gInfo);
		w2 = new RealKernel(gInfo);
		w3 = new RealKernel(gInfo);
		w1v = new RealKernel(gInfo);
		w2m = new RealKernel(gInfo);
		ErfFMTweight erfFMTweight(sphereRadius, sphereSigma);
		applyFuncGsq(gInfo, erfFMTweight, w0->data, w1->data, w2->data, w3->data, w1v->data, w2m->data);
		w0->set();
		w1->set();
		w2->set();
		w3->set();
		w1v->set();
		w2m->set();
	}
}

SiteProperties::~SiteProperties()
{	if(sphereRadius)
	{	delete w0;
		delete w1;
		delete w2;
		delete w3;
		delete w1v;
		delete w2m;
	}
}

double Molecule::get_charge() const
{	double Q = 0.0;
	for(int i=0; i<nSites; i++)
	{	SiteProperties& s = *site[i].prop;
		if(s.chargeZ && s.chargeKernel)
			Q += s.chargeZ * s.chargeKernel->data[0];
	}
	return Q;
}

double Molecule::get_dipole() const
{	vector3<> electricP(0,0,0);
	for(int i=0; i<nSites; i++)
	{	SiteProperties& s = *site[i].prop;
		if(s.chargeZ && s.chargeKernel)
			electricP += site[i].pos * s.chargeZ * s.chargeKernel->data[0];
	}
	//Check that dipole (if any) is lined up with z-axis
	double dipoleOffZaxis = hypot(electricP[0], electricP[1]);
	if(dipoleOffZaxis > std::max(1e-6, 1e-10*electricP.length()))
		die("Fluid molecule dipole moment component off z-axis is %lg.\n Please orient fluid dipole along reference z-axis\n", dipoleOffZaxis);
	return electricP[2];
}

std::map<double,int> Molecule::getBonds() const
{	std::map<double,int> bond;
	for(int i=0; i<nSites; i++)
	{	double Ri = site[i].prop->sphereRadius;
		if(Ri)
		{	for(int j=i+1; j<nSites; j++)
			{	double Rj = site[j].prop->sphereRadius;
				if(Rj)
				{	if(fabs(Ri+Rj-(site[i].pos-site[j].pos).length()) < 1e-6*(Ri+Rj))
					{	//In contact within tolerance:
						bond[Ri*Rj/(Ri+Rj)]++;
					}
				}
			}
		}
	}
	return bond;
}
