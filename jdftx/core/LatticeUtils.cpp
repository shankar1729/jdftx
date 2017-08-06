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

#include <core/LatticeUtils.h>
#include <core/WignerSeitz.h>
#include <core/Util.h>
#include <core/matrix.h>
#include <algorithm>
#include <cfloat>
#include <list>

double symmThreshold = 1e-4;
double symmThresholdSq = symmThreshold*symmThreshold;

//Compute reduced lattice vectors:
matrix3<> reduceLatticeVectors(const matrix3<>& R, matrix3<int>* transmission, matrix3<int>* invTransmission)
{	matrix3<> Rreduced = R;
	matrix3<int> t (1,1,1);
	matrix3<int> tInv(1,1,1);
	while(true)
	{	bool changed = false;
		for(int k1 = 0; k1 < 3; k1 ++)
		{	int k2 = (k1+1)%3;
			int k3 = (k1+2)%3;
			for(int i=-1; i<=1; i++)
				for(int j=-1; j<=1; j++)
				{	//Add/subtract upto one each of k2 and k3'th directions to the k1st direction:
					matrix3<int> d(1,1,1), dInv(1,1,1);
					d(k2,k1)=i; d(k3,k1)=j; 
					dInv(k2,k1)=-i; dInv(k3,k1)=-j; 
					
					//Check if that transformation reduces R:
					matrix3<> Rproposed = Rreduced * d;
					if(nrm2(Rproposed) < nrm2(Rreduced) - symmThreshold)
					{	changed = true;
						Rreduced = Rproposed;
						t = t * d;
						tInv = dInv * tInv;
					}
				}
		}
		if(!changed) break;
	}
	if(transmission) *transmission = t;
	if(invTransmission) *invTransmission = tInv;
	return Rreduced;
}

//Determine symmetries of Bravais lattice:
std::vector<matrix3<int>> getSymmetries(const matrix3<>& R, vector3<bool> isTruncated,
	matrix3<>* Rreduced, matrix3<int>* transmission, matrix3<int>* invTransmission)
{
	//Get reduced lattice vectors:
	matrix3<int> t, tInv;
	matrix3<> Rred = reduceLatticeVectors(R, &t, &tInv);
	
	//Check symmetries by finding integer matrices that leave the metric invariant:
	matrix3<> metric = (~Rred)*Rred;
	std::vector<matrix3<int>> sym;
	//loop over integer components for each matrix entry:
	//all entries of reduced lattice vector symmetries must be -1, 0 or 1.
	matrix3<int> m;
	#define iLOOP(x) for(x=-1; x<=1; x++)
	iLOOP(m(0,0)) iLOOP(m(0,1)) iLOOP(m(0,2))
	iLOOP(m(1,0)) iLOOP(m(1,1)) iLOOP(m(1,2))
	iLOOP(m(2,0)) iLOOP(m(2,1)) iLOOP(m(2,2))
		if(nrm2(metric - (~m)*metric*m) < symmThreshold * nrm2(metric)) //leaves metric invariant
		{	matrix3<int> mOrig = t * m * tInv; //in original lattice coordinates
			//Check if symmetry connects truncated and untruncated directions:
			bool symBroken = false;
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
					if(mOrig(i,j) && isTruncated[i]!=isTruncated[j])
						symBroken = true;
			if(!symBroken) sym.push_back(mOrig);
		}
	#undef  iLOOP
	
	if(Rreduced) *Rreduced = Rred;
	if(transmission) *transmission = t;
	if(invTransmission) *invTransmission = tInv;
	return sym;
}

Supercell::Supercell(const GridInfo& gInfo,
	const std::vector<vector3<>>& kmeshReduced,
	const std::vector<SpaceGroupOp>& sym, const std::vector<int>& invertList)
: gInfo(gInfo)
{
	logPrintf("\n----- Initializing Supercell corresponding to k-point mesh -----\n");

	//Compute kmesh = closure of kmeshReduced under symmetry group, sym:
	PeriodicLookup< vector3<> > plook(kmesh, gInfo.GGT, kmeshReduced.size()*sym.size()); //look-up table for O(1) fuzzy searching
	for(int invert: invertList)
		for(unsigned iReduced=0; iReduced<kmeshReduced.size(); iReduced++)
		{	const vector3<>& kOrig = kmeshReduced[iReduced];
			for(unsigned iSym=0; iSym<sym.size(); iSym++)
			{	const matrix3<int>& m = sym[iSym].rot;
				vector3<> k = (~m) * kOrig * invert;
				//Reduce to centered zone (each reciprocal lattice coord in [-0.5,0.5))
				vector3<int> offset;
				for(int i=0; i<3; i++)
				{	offset[i] = -floor(0.5+k[i]);
					k[i] += offset[i];
				}
				//Add to map if this k-vector has not yet been encountered:
				if(plook.find(k) == string::npos)
				{	plook.addPoint(kmesh.size(), k);
					kmesh.push_back(k);
					KmeshTransform kTransform = { iReduced, iSym, invert, offset };
					kmeshTransform.push_back(kTransform);
				}
			}
		}
	
	//Construct a linearly-independent basis for the k-points:
	matrix3<> kBasis;
	{	//Add periodic images from neighbouring Brillouin zones:
		std::vector<vector3<>> kmesh333 = kmesh;
		vector3<int> iZone;
		for(iZone[0]=-1; iZone[0]<=1; iZone[0]++)
		for(iZone[1]=-1; iZone[1]<=1; iZone[1]++)
		for(iZone[2]=-1; iZone[2]<=1; iZone[2]++)
			if(iZone.length_squared())
				for(vector3<> k: kmesh)
					kmesh333.push_back(k + iZone);
		vector3<> v[3];
		//Pick v[0] to be the shortest dk:
		double minLengthSq = DBL_MAX;
		for(size_t i=1; i<kmesh333.size(); i++)
		{	vector3<> dk = kmesh333[i] - kmesh333.front();
			double lengthSq = dk.length_squared();
			if(lengthSq < minLengthSq)
			{	v[0] = dk;
				minLengthSq = lengthSq;
			}
		}
		//Pick v[1] to be the shortest dk not linearly dependent with v[0]:
		minLengthSq = DBL_MAX;
		double crossThreshold = symmThreshold * v[0].length();
		for(size_t i=1; i<kmesh333.size(); i++)
		{	vector3<> dk = kmesh333[i] - kmesh333.front();
			if(cross(v[0], dk).length() > crossThreshold)
			{	double lengthSq = dk.length_squared();
				if(lengthSq < minLengthSq)
				{	v[1] = dk;
					minLengthSq = lengthSq;
				}
			}
		}
		//Pick v[1] to be the shortest dk not linearly dependent with v[0] and v[1]:
		minLengthSq = DBL_MAX;
		double boxThreshold = symmThreshold * cross(v[0],v[1]).length();
		for(size_t i=1; i<kmesh333.size(); i++)
		{	vector3<> dk = kmesh333[i] - kmesh333.front();
			if(fabs(box(v[0], v[1], dk)) > boxThreshold)
			{	double lengthSq = dk.length_squared();
				if(lengthSq < minLengthSq)
				{	v[2] = dk;
					minLengthSq = lengthSq;
				}
			}
		}
		for(int k=0; k<3; k++) kBasis.set_col(k, v[k]);
	}
	
	//Check if kmesh forms a Bravais lattice:
	kBasis = reduceLatticeVectors(kBasis);
	matrix3<> kBasisInv = inv(kBasis);
	//--- check that each kpoint is integral in above basis
	#define kMeshErrorMsg \
		"k-point mesh is not a Bravais lattice, which is currently required for\n" \
		"exact-exchange and density-of-states calculations.\n" \
		"HINT: This might be caused by an off-Gamma k-point mesh (eg. Monkhorst-Pack)\n" \
		"   whose symmetries are a subgroup of that of the system. Try switching to\n" \
		"   Gamma-centered sampling or disabling symmetries.\n"
	for(vector3<> kpoint: kmesh)
	{	double err; round(kBasisInv * (kpoint-kmesh.front()), &err);
		if(err > symmThreshold) die(kMeshErrorMsg);
	}
	//--- check that the k-mesh covers the Brillouin zone
	if(fabs(fabs(det(kBasisInv)) - kmesh.size()) > kmesh.size() * symmThreshold)
		die(kMeshErrorMsg);
	#undef kMeshErrorMsg
	
	//Compute the supercell matrix:
	matrix3<> superTemp = inv(gInfo.R) * reduceLatticeVectors(gInfo.R * ~kBasisInv);
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{	super(i,j) = int(round(superTemp(i,j)));
			if(fabs(superTemp(i,j) - super(i,j)) > symmThreshold)
				die("k-point mesh does not correspond to a commensurate (integer)\n"
					"super-cell. This is currently required for exact-exchange\n"
					"and density-of-states calculations.\n");
		}
	
	//Pivot columns to get closest to identity:
	{	int permutation[6][3] = { {0,1,2}, {1,2,0}, {2,0,1}, {2,1,0}, {1,0,2}, {0,2,1} };
		double maxScore = 0; int pBest = 0;
		for(int p=0; p<6; p++)
		{	double score = 1.;
			for(int k=0; k<3; k++)
			{	vector3<int> col = super.column(permutation[p][k]);
				score *= pow(col[k], 2);
				score /= col.length_squared();
			}
			if(score > maxScore)
			{	maxScore = score;
				pBest = p;
			}
		}
		matrix3<int> superOrig = super;
		for(int k=0; k<3; k++)
			super.set_col(k, superOrig.column(permutation[pBest][k]));
	}
	if(det(super) < 0) super *= -1;
	logPrintf("Lattice vector linear combinations in supercell:\n");
	super.print(globalLog, " %2d ");
	Rsuper = gInfo.R * super;
	logPrintf("Supercell lattice vectors:\n");
	Rsuper.print(globalLog, " %lg ");
}


std::map<vector3<int>, matrix> getCellMap(const matrix3<>& R, const matrix3<>& Rsup, const vector3<bool>& isTruncated,
	const std::vector<vector3<>>& x1, const std::vector<vector3<>>& x2, double rSmooth, string fname)
{
	//Get neighbouring supercell lattice points defining supercell WS cell:
	logSuspend();
	WignerSeitz wsSup(Rsup);
	logResume();
	std::vector<vector3<int>> xSupNeigh = wsSup.getNeighbours(isTruncated);
	std::vector<vector3<>> rSupNeigh; //Cartesian version of above
	for(const vector3<int>& xSup: xSupNeigh)
		rSupNeigh.push_back(Rsup * xSup);
	
	//Calculate range of neighbouring cells required:
	double d12sqMax = 0.;
	matrix3<> RTR = (~R)*R;
	for(const vector3<>& x1i: x1)
		for(const vector3<>& x2j: x2)
		{	vector3<> dx = x1i - x2j;
			for(int l=0; l<3; l++) if(isTruncated[l]) dx[l] = 0.; //project out truncated directions
			d12sqMax = std::max(d12sqMax, RTR.metric_length_squared(dx));
		}
	double Rmax = wsSup.circumRadius() + sqrt(d12sqMax);
	
	//Collect lattice vectors and weight matrices:
	matrix3<> invR = inv(R);
	vector3<int> iCellMax; //bounding box of rMax (in unit-cell lattice coordinates)
	for(int l=0; l<3; l++)
		iCellMax[l] = isTruncated[l] ? 0 : (1 + int(ceil(Rmax * invR.row(l).length())));
	vector3<int> iCell;
	std::map<vector3<int>, matrix> iCellMap;
	for(iCell[0]=-iCellMax[0]; iCell[0]<=iCellMax[0]; iCell[0]++)
	for(iCell[1]=-iCellMax[1]; iCell[1]<=iCellMax[1]; iCell[1]++)
	for(iCell[2]=-iCellMax[2]; iCell[2]<=iCellMax[2]; iCell[2]++)
	{	matrix w = zeroes(x1.size(), x2.size());
		bool nonZero = false;
		for(size_t i=0; i<x1.size(); i++)
		for(size_t j=0; j<x2.size(); j++)
		{	vector3<> dx = iCell + x1[i] -x2[j];
			for(int l=0; l<3; l++) if(isTruncated[l]) dx[l] = 0.; //project out truncated directions
			vector3<> dr = R*dx; //cartesian displacement
			double drSq0 = dr.length_squared(); //distance from origin
			//Get distance to nearest supercell neighbour:
			double drSqNeigh = DBL_MAX;
			for(const vector3<>& rSup: rSupNeigh)
				drSqNeigh = std::min(drSqNeigh, (dr-rSup).length_squared());
			//Calculate weight:
			double t = 0.5*(drSq0-drSqNeigh)/rSmooth; //<0 inside and >0 outside
			if(t < 1.) //within smoothed margin of wsSup
			{	w.set(i,j, (t<=-1.) ? 1. : 0.25*(2.-t*(3.-t*t)));
				nonZero = true;
			}
		}
		if(nonZero) iCellMap[iCell] = w;
	}
	
	//Normalize weights:
	std::set<vector3<int>> cellsDone;
	typedef std::map<vector3<int>, matrix>::iterator Iter;
	matrix3<> superInv = inv(Rsup) * R;
	for(Iter iter=iCellMap.begin(); iter!=iCellMap.end(); iter++) if(!cellsDone.count(iter->first))
	{	std::vector<Iter> equiv(1, iter);
		matrix wSum = iter->second;
		//Loop over equivalent cells:
		Iter iter2=iter; iter2++;
		for(; iter2!=iCellMap.end(); iter2++) if(!cellsDone.count(iter2->first))
		{	double err; round(superInv * (iter2->first - iter->first), &err); //difference in super-lattice coordinates
			if(err < symmThreshold) //integer superlattice offset
			{	equiv.push_back(iter2);
				wSum += iter2->second;
				cellsDone.insert(iter2->first); //mark as done
			}
		}
		//Replace weight sum with reciprocal:
		for(size_t i=0; i<x1.size(); i++)
		for(size_t j=0; j<x2.size(); j++)
		{	double& wCur = wSum.data()[wSum.index(i,j)].real();
			assert(wCur > symmThreshold);
			wCur = 1./wCur;
		}
		//Normalize weights element-wise:
		for(Iter& iter2: equiv) scale(wSum, iter2->second);
	}

	//Write the cell map if requested
	if(mpiUtil->isHead() && fname.length())
	{	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "#i0 i1 i2  x y z  (integer lattice combinations, and cartesian offsets)\n");
		for(const auto& entry: iCellMap)
		{	const vector3<int>& i = entry.first;
			vector3<> r = R * i;
			fprintf(fp, "%+2d %+2d %+2d  %+11.6lf %+11.6lf %+11.6lf\n", i[0], i[1], i[2], r[0], r[1], r[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	return iCellMap;
}
