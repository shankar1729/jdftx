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

#include <core/WignerSeitz.h>
#include <core/Util.h>

static const double minDistSq = 1e-16; //threshold on distance squared
static const double geomRelTol = 1e-14; //relative tolerance to geometry (tesselation volumes, orthogonality checks etc)

//Construct Wigner-Seitz cell given lattice vectors
WignerSeitz::WignerSeitz(const matrix3<>& R) : R(R), invR(inv(R)), RTR((~R)*R)
{	logPrintf("Constructing Wigner-Seitz cell: "); logFlush();
	//Initialize with parallelopiped formed from nearest neighbours:
	matrix3<> eqns = matrix3<>(2./RTR(0,0), 2./RTR(1,1), 2./RTR(2,2)) * RTR;
	matrix3<> eqnsInv = inv(eqns);
	Vertex* vInit[2][2][2];
	double rMaxSq = 0.; //current maximum distance squared of cell edges from origin
	for(int i0=0; i0<2; i0++)
		for(int i1=0; i1<2; i1++)
			for(int i2=0; i2<2; i2++)
			{	Vertex* vNew = new Vertex;
				vInit[i0][i1][i2] = vNew;
				vertex.push_back(vNew);
				vNew->pos = eqnsInv * vector3<int>(2*i0-1, 2*i1-1, 2*i2-1);
				double rSq = RTR.metric_length_squared(vNew->pos);
				if(rSq > rMaxSq) rMaxSq = rSq;
			}
	for(int dir=0; dir<3; dir++) //normal direction
		for(int s=0; s<2; s++) //sign for direction
		{	Face* fNew = new Face;
			face.insert(fNew);
			vector3<int> i;
			int xDir = s ? (dir+1)%3 : (dir+2)%3;
			int yDir = s ? (dir+2)%3 : (dir+1)%3;
			i[dir] = 2*s-1;
			i[xDir] = 0; i[yDir] = 0;
			fNew->img = i;
			fNew->eqn = (2*s-1) * eqns.row(dir);
			i[dir] = s;
			i[xDir] = 0; i[yDir] = 0; Vertex* vA = vInit[i[0]][i[1]][i[2]];
			i[xDir] = 1; i[yDir] = 0; Vertex* vB = vInit[i[0]][i[1]][i[2]];
			i[xDir] = 1; i[yDir] = 1; Vertex* vC = vInit[i[0]][i[1]][i[2]];
			i[xDir] = 0; i[yDir] = 1; Vertex* vD = vInit[i[0]][i[1]][i[2]];
			addEdge(fNew, vA, vB);
			addEdge(fNew, vB, vC);
			addEdge(fNew, vC, vD);
			addEdge(fNew, vD, vA);
		}
	
	//Determine coordinate range to explore:
	vector3<int> ivMax;
	matrix3<> invR = inv(R);
	for(int dir=0; dir<3; dir++)
		ivMax[dir] = ceil(2. * sqrt(rMaxSq * invR.row(dir).length_squared()));
	//List lattice points to try in increasing order of distance:
	std::multimap<double, vector3<int> > ivMap;
	vector3<int> iv;
	for(iv[0]=-ivMax[0]; iv[0]<=ivMax[0]; iv[0]++)
	for(iv[1]=-ivMax[1]; iv[1]<=ivMax[1]; iv[1]++)
	for(iv[2]=-ivMax[2]; iv[2]<=ivMax[2]; iv[2]++)
		if(iv.length_squared() > 1) //exclude the planes already considered
		{	double rSq = RTR.metric_length_squared(iv);
			if(rSq < 4*rMaxSq)
				ivMap.insert(std::make_pair(rSq, iv));
		}
	
	//Slice the current zone in increasing order of point distance:
	for(const std::pair<double, vector3<int> >& ivMapEntry: ivMap)
	{	if(ivMapEntry.first > 4*rMaxSq) break; //in that case bisector cannot touch surface
		//Slice with plane:
		addPlane(ivMapEntry.second);
		//Update rMaxSq (cell gets smaller as more slice planes are added)
		rMaxSq = 0.;
		for(const Vertex* v: vertex)
		{	double rSq = RTR.metric_length_squared(v->pos);
			if(rSq > rMaxSq) rMaxSq = rSq;
		}
	}
	
	//Pick one face out of each set of inversion partners:
	std::map<vector3<int>,Face*> imgFace; //map from Face::img to Face*
	for(Face* f: face)
	{	//Search for inversion partner:
		auto imgIter = imgFace.find(-f->img);
		if(imgIter == imgFace.end()) //no previous inversion partner
		{	imgFace[f->img] = f; //add to inversion partner map
			faceHalf.push_back(f); //add face to half face list
		}
		else //found inversion partner
		{	//Check inversion symmetry
			Face* fPair[2] = {f, imgIter->second};
			std::list<vector3<>> vPair[2];
			//Extract vertices for both faces:
			for(unsigned i=0; i<2; i++)
				for(const Edge* e: fPair[i]->edge)
					vPair[i].push_back(e->vertex[e->face[0]==fPair[i] ? 0 : 1]->pos);
			assert(vPair[0].size() == vPair[1].size());
			//Align the vertex lists by finding inversion partner for first vertex of first list
			vector3<> v0start = vPair[0].front();
			auto v1iter = vPair[1].begin();
			for(; v1iter != vPair[1].end(); v1iter++)
				if(RTR.metric_length_squared(*v1iter + v0start) < minDistSq)
					break;
			assert(v1iter != vPair[1].end());
			//Check symmetry of each pair of vertices:
			auto v0iter = vPair[0].begin();
			while(v0iter != vPair[0].end())
			{	assert(RTR.metric_length_squared(*v0iter + *v1iter) < minDistSq);
				if(v1iter == vPair[1].begin()) v1iter = vPair[1].end();
				v0iter++; v1iter--; //note vertices in opposite order due to inversion
			}
			imgFace.erase(imgIter); //check off pair from inversion partner map
		}
	}
	assert(imgFace.size() == 0); //Check that all faces occur in inversion symmetric pairs
	
	//Print summary:
	int nQuad=0, nHex=0;
	for(const Face* f: faceHalf)
	{	if(f->edge.size()==4) nQuad += 2;
		if(f->edge.size()==6) nHex += 2;
	}
	logPrintf("%lu faces (%d quadrilaterals, %d hexagons)\n", face.size(), nQuad, nHex);
}

WignerSeitz::~WignerSeitz()
{	for(Vertex* v: vertex) delete v;
	for(Edge* e: edge) delete e;
	for(Face* f: face) delete f;
}

//Output a list of simplexes that tesselate half the Weigner-Seitz cell (remainder filled by inversion symmetry)
std::vector<Simplex<3>> WignerSeitz::getSimplices() const
{	const double Vcell = fabs(det(R)); //expected volume of Wigner-Seitz cell
	const double Vtol = geomRelTol * Vcell;
	std::vector<Simplex<3>> sArr;
	double Vtot = 0.;
	for(const Face* f: faceHalf)
	{	//Collect vertices in order:
		std::vector<Simplex<3>::Point> v;
		for(const Edge* e: f->edge)
		{	vector3<> vCart = R * e->vertex[(e->face[0]==f) ? 0 : 1]->pos; //in cartesian coordinates
			v.push_back({{vCart[0], vCart[1], vCart[2]}}); //convert to Simplex<3>::Point
		}
		//Collect simplices from this face:
		for(unsigned i=2; i<v.size(); i++)
		{	Simplex<3> s;
			s.v[0] = v[0];
			s.v[1] = v[i-1];
			s.v[2] = v[i];
			s.init();
			Vtot += (2./6) * s.V; //Factor of 2 includes inversion partner
			sArr.push_back(s);
		}
	}
	assert(fabs(Vtot-Vcell) < Vtol); //Check total volume
	return sArr;
}

//Output a list of simplexes that tesselate half the 2D Weigner-Seitz cell
//(remainder filled by inversion symmetry) after ignoring lattice direction iDir
std::vector<Simplex<2>> WignerSeitz::getSimplices(int iDir) const
{	matrix3<> Rplanar = getRplanar(iDir);
	const double Acell = fabs(det(R)) / R.column(iDir).length(); //expected area of 2D Wigner-Seitz cell
	const double Atol = geomRelTol * Acell;
	std::vector<Simplex<2>> sArr;
	double Atot = 0.;
	for(const Face* f: faceHalf)
		if(f->img[iDir]==1 && f->img.length_squared()==1)
		{	//Collect vertices in order:
			std::vector<Simplex<2>::Point> v;
			for(const Edge* e: f->edge)
			{	vector3<> vCart = Rplanar * e->vertex[(e->face[0]==f) ? 0 : 1]->pos; //in rotated cartesian coordinates
				v.push_back({{vCart[0], vCart[1]}}); //convert to Simplex<2>::Point, dropping iDir (which is z after rotation)
			}
			//Add simplices:
			for(unsigned i=0; i<v.size(); i++)
			{	Simplex<2> s;
				s.v[0] = v[i ? i-1 : v.size()-1];
				s.v[1] = v[i];
				s.init();
				Atot += s.V;
				sArr.push_back(s);
			}
			break;
		}
	assert(fabs(Atot-Acell) < Atol); //Check total volume
	return sArr;
}

//Rotated lattice vectors such that iDir maps onto z-axis
matrix3<> WignerSeitz::getRplanar(int iDir) const
{	int jDir = (iDir+1)%3;
	int kDir = (iDir+2)%3;
	//Construct an orthogonal matrix that takes z-axis to iDir and x-axis to jDir
	matrix3<> rot;
	rot.set_col(2, R.column(iDir) / R.column(iDir).length());
	rot.set_col(0, R.column(jDir) / R.column(jDir).length());
	assert(fabs(dot(rot.column(2), rot.column(0))) < geomRelTol);
	assert(fabs(dot(rot.column(2), R.column(kDir))) < geomRelTol * R.column(kDir).length());
	rot.set_col(1, cross(rot.column(2), rot.column(0)));
	//Apply inverse rotation to R to get new lattice vectors:
	return (~rot) * R;
}


//Radius of largest inscribable sphere (circle) centered at origin
double WignerSeitz::inRadius(int iDir) const
{	double rSqMin = DBL_MAX;
	for(const Face* f: faceHalf)
		if(iDir<0 || !f->img[iDir]) //ignore direction iDir, if non-negative
		{	double rSq = RTR.metric_length_squared(f->img);
			if(rSq < rSqMin) rSqMin = rSq;
		}
	return 0.5*sqrt(rSqMin);
}

//Radius of smallest circumscribable sphere (circle) centered at origin
double WignerSeitz::circumRadius(int iDir) const
{	double rSqMax = 0.;
	for(const Vertex* v: vertex)
	{	vector3<> pos = v->pos;
		if(iDir>=0) pos[iDir] = 0.; //project out iDir, if non-negative
		double rSq = RTR.metric_length_squared(pos);
		if(rSq > rSqMax) rSqMax = rSq;
	}
	return sqrt(rSqMax);
}

//Write a wireframe plot to file (for gnuplot)
void WignerSeitz::writeWireframePlot(const char* filename) const
{	FILE* fp = fopen(filename, "w");
	for(const Face* f: face)
	{	bool polyStarted = false;
		for(const Edge* e: f->edge)
		{	int iStart = (e->face[0]==f) ? 0 : 1;
			vector3<> startPos = R * e->vertex[iStart]->pos;
			vector3<> endPos = R * e->vertex[1-iStart]->pos;
			if(!polyStarted)
			{	polyStarted = true;
				fprintf(fp, "%lf\t%lf\t%lf\n", startPos[0], startPos[1], startPos[2]);
			}
			fprintf(fp, "%lf\t%lf\t%lf\n", endPos[0], endPos[1], endPos[2]);
		}
		fprintf(fp, "\n\n");
	}
	fclose(fp);
}

//Write a wireframe plot for DataExplorer
void WignerSeitz::writeWireframeDX(const char* filename) const
{	FILE* fp = fopen(filename, "w");
	//Output vertices:
	fprintf(fp, "object \"pos\" array items %lu\n", vertex.size());
	fprintf(fp, "\trank 1 type float shape 3\n");
	fprintf(fp, "\tdata follows\n");
	int vid = 0; std::map<const Vertex*, int> vID; //vertex IDs
	for(const Vertex* v: vertex)
	{	vector3<> pos = R * v->pos; //cartesian coords
		fprintf(fp, "\t%lg\t%lg\t%lg\n", pos[0], pos[1], pos[2]);
		vID[v] = vid++; //remember index for each pointer
	}
	fprintf(fp, "\n");
	//Output edges:
	fprintf(fp, "object \"conn\" class array type int rank 1 shape 2 items %lu\n", edge.size());
	fprintf(fp, "\tattribute \"element type\" string \"lines\"\n");
	fprintf(fp, "\tattribute \"ref\" string \"positions\"\n");
	fprintf(fp, "\tdata follows\n");
	for(const Edge* e: edge)
		fprintf(fp, "\t%d\t%d\n", vID[e->vertex[0]], vID[e->vertex[1]]);
	fprintf(fp, "\n");
	//Put together wire mesh object:
	fprintf(fp, "object \"atoms\" field\n");
	fprintf(fp, "\tcomponent \"positions\" \"pos\"\n");
	fprintf(fp, "\tcomponent \"connections\" \"conn\"\n");
	fclose(fp);
}


//Check if the data structure is valid (all links are reciprocated etc.)
void WignerSeitz::checkGraph() const
{	for(const Vertex* v: vertex)
		for(const Edge* e: v->edge)
			if(e->vertex[0] != v && e->vertex[1] != v)
				die("Vertex (%lf,%lf,%lf) links to an edge that does not link back to it.\n", v->pos[0], v->pos[1], v->pos[2])
	for(const Edge* e: edge)
	{	for(const Vertex* v: e->vertex)
			if(std::find(v->edge.begin(), v->edge.end(), e) == v->edge.end())
				die("Vertex (%lf,%lf,%lf) does not link to an edge that links to it.\n", v->pos[0], v->pos[1], v->pos[2])
		for(const Face* f: e->face)
			if(std::find(f->edge.begin(), f->edge.end(), e) == f->edge.end())
				die("Face (%d,%d,%d) does not link to an edge that links to it.\n", f->img[0], f->img[1], f->img[2])
	}
	for(const Face* f: face)
		for(const Edge* e: f->edge)
			if(e->face[0] != f && e->face[1] != f)
				die("Face (%d,%d,%d) links to an edge that does not link back to it.\n", f->img[0], f->img[1], f->img[2])
}

//Output vertex, edge and face connectivity info:
void WignerSeitz::writeGraph(FILE* fp) const
{	//Output vertices:
	std::map<const Vertex*,int> vID;
	fprintf(fp, "--- Vertices ---\n");
	int vid = 0;
	for(const Vertex* v: vertex)
		fprintf(fp, "%2d: (%lf, %lf, %lf)\n", (vID[v] = ++vid), v->pos[0], v->pos[1], v->pos[2]);
	//Output faces as list of edges with vertex numbers:
	fprintf(fp, "--- Faces ---\n");
	for(const Face* f: face)
	{	for(const Edge* e: f->edge)
		{	int iStart = e->face[0]==f ? 0 : 1;
			fprintf(fp, "(%d->%d) ", vID[e->vertex[iStart]], vID[e->vertex[1-iStart]]);
		}
		fprintf(fp, "\n");
	}
}

//Slice the current polyhedron by the perpendicular bisector of 0->a (a in lattice coordinates)
void WignerSeitz::addPlane(const vector3<int>& a)
{	vector3<> eqn = (2./RTR.metric_length_squared(a))*(RTR*a); //eqn.x==1 is perpendicular bisector plane
	//Check each vertex for removal:
	std::set<Vertex*> onPlane; //set of vertices on new plane
	for(std::list<Vertex*>::iterator vIter=vertex.begin(); vIter!=vertex.end();)
	{	Vertex* v = *vIter;
		double d = dot(eqn, v->pos) - 1.; //projected distance of vertex (0 on plane, <0 interior)
		if(d > 0.) //vertex on the exterior => remove:
		{	//Trim each edge ending on this vertex:
			for(Edge* e: v->edge)
			{	//Check other vertex on this edge
				int iOther = (e->vertex[0]==v) ? 1 : 0;
				double dOther = dot(eqn, e->vertex[iOther]->pos) - 1.;
				if(dOther > 0.) //remove edge:
				{	e->vertex[iOther]->edge.remove(e); //remove from other vertex
					for(Face* f: e->face) f->edge.remove(e); //remove from faces
					edge.erase(e); //remove from global list
					delete e;
				}
				else //trim edge:
				{	Vertex* vNew = new Vertex;
					vNew->pos = v->pos + (d/(d-dOther))*(e->vertex[iOther]->pos - v->pos);
					vNew->edge.push_back(e); //add e to edges ending on new vertex
					e->vertex[1-iOther] = vNew; //replace v by vNew in edge e
					onPlane.insert(vNew); //add to list of vertices on plane
				}
			}
			//vertex unlinked from edges, can remove it now:
			vIter = vertex.erase(vIter); //remove from global list
			delete v;
		}
		else vIter++; //vertex on boundary or interior => retain:
	}
	vertex.insert(vertex.end(), onPlane.begin(), onPlane.end()); //add new vertices to global list
	//Create a map from sliced faces to their vertices on the slice plane:
	std::map<Face*, std::set<Vertex*> > sliceVert; 
	for(Vertex* v: onPlane)
		for(Edge* e: v->edge)
			for(Face* f: e->face)
				sliceVert[f].insert(v);
	for(auto i=sliceVert.begin(); i!=sliceVert.end();)
	{	switch(i->second.size())
		{	case 1: sliceVert.erase(i++); break; //remove single vertex intersections and proceed
			case 2: i++; break; //case of interest: intersection at two points (proceed)
			default: assert(!"Error: face sliced at more than 2 points."); //should never happen
		}
	}
	assert(sliceVert.size() == onPlane.size());
	if(!sliceVert.size()) return;
	//Complete sliced faces in sequence (based on shared vertices), and create new face:
	Face* fNew = new Face;
	fNew->img = a;
	fNew->eqn = eqn;
	face.insert(fNew);
	Vertex* vPrev = 0;
	auto curSlice = sliceVert.begin();
	while(sliceVert.size())
	{	Face* f = curSlice->first;
		if(!vPrev) //Pick the starting vertex after deciding polygon direction
		{	//Get the two vertices for this face slice: 
			Vertex* vA = *curSlice->second.begin();
			Vertex* vB = *curSlice->second.rbegin();
			//Determine order (note order in new face is opposite from that in old):
			for(Edge* e: f->edge)
			{	int iEnd = (e->face[0]==f) ? 1 : 0;
				if(e->vertex[iEnd]==vA) { vPrev=vB; break; }
				if(e->vertex[iEnd]==vB) { vPrev=vA; break; }
			}
		}
		curSlice->second.erase(vPrev);
		Vertex* v = *curSlice->second.begin();
		addEdge(f, v, vPrev, true); //complete face f (note opposite order)
		addEdge(fNew, vPrev, v);
		sliceVert.erase(curSlice);
		//Search for other old face that conatins v:
		bool nextSliceFound = false;
		for(auto nextSlice=sliceVert.begin(); nextSlice!=sliceVert.end(); nextSlice++)
			if(nextSlice->second.count(v))
			{	curSlice = nextSlice;
				vPrev = v;
				nextSliceFound = true;
			}
		assert(!sliceVert.size() || nextSliceFound);
	}
	
	//Weld almost identical vertices (should have been identical in exact arithmetic): [version 2]
	for(auto vIter=vertex.begin(); vIter!=vertex.end(); vIter++)
	{	Vertex* v = *vIter;
		for(auto eIter=v->edge.begin(); eIter!=v->edge.end();)
		{	Edge* e = *eIter;
			Vertex* vOther = e->vertex[(e->vertex[0]==v) ? 1 : 0];
			if(RTR.metric_length_squared(v->pos - vOther->pos) > minDistSq)
				eIter++; //do nothing
			else //Found a neighbour vOther within tolerance:
			{	//Remove the trivial edge between v and vOther:
				eIter = v->edge.erase(eIter); //remove from v
				vOther->edge.remove(e); //remove from vOther
				for(Face* f: e->face) f->edge.remove(e); //remove from faces
				edge.erase(e); //remove from global list
				delete e;
				//Replace vOther with v, and get rid of vOther:
				for(Edge* eOther: vOther->edge)
				{	eOther->vertex[(eOther->vertex[0]==vOther) ? 0 : 1] = v;
					v->edge.push_back(eOther);
				}
				vertex.remove(vOther);
				delete vOther;
			}
		}
	}
	//Remove null faces (could have been caused by welding)
	for(auto fIter=face.begin(); fIter!=face.end();)
		if(!(*fIter)->edge.size())
			fIter = face.erase(fIter);
		else fIter++;

}

//Add a edge from vStart towards vEnd in face f
void WignerSeitz::addEdge(Face* f, Vertex* vStart, Vertex* vEnd, bool lastEdge)
{	Edge* e = 0;
	for(Edge* eOld: edge)
		if(eOld->vertex[1]==vStart && eOld->vertex[0]==vEnd)
		{	e = eOld;
			break;
		}
		else assert(!(eOld->vertex[0]==vStart && eOld->vertex[1]==vEnd)); //existing edge should not have same direction
	if(e) //Link existing edge to face:
	{	assert(!e->face[1]); //edge should not be associated with two faces already
		e->face[1] = f;
	}
	else //Create new edge:
	{	e = new Edge;
		edge.insert(e);
		e->vertex[0] = vStart;
		e->vertex[1] = vEnd;
		e->face[0] = f;
		e->face[1] = 0;
		vStart->edge.push_back(e);
		vEnd->edge.push_back(e);
	}
	//Add edge to face:
	if(lastEdge) //Find insert location:
	{	auto eIter = f->edge.begin();
		while(eIter != f->edge.end())
		{	Edge* ePrev = *(eIter++);
			Vertex* vPrev = (ePrev->face[0]==f) ? ePrev->vertex[1] : ePrev->vertex[0]; //select ending vertex by direction
			if(vPrev == vStart) break; //eIter now points to the location to insert new edge at
		}
		f->edge.insert(eIter, e);
	}
	else f->edge.push_back(e);
}
