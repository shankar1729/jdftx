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

#include <map>
#include <string>

using namespace std;

const int nElements=118;
pair<string,int> atomicSymbolArray[nElements] =
{
	{"Ac", 89},
	{"Ag", 47},
	{"Al", 13},
	{"Am", 95},
	{"Ar", 18},
	{"As", 33},
	{"At", 85},
	{"Au", 79},
	{"B", 5},
	{"Ba", 56},
	{"Be", 4},
	{"Bh", 107},
	{"Bi", 83},
	{"Bk", 97},
	{"Br", 35},
	{"C", 6},
	{"Ca", 20},
	{"Cd", 48},
	{"Ce", 58},
	{"Cf", 98},
	{"Cl", 17},
	{"Cm", 96},
	{"Cn", 112},
	{"Co", 27},
	{"Cr", 24},
	{"Cs", 55},
	{"Cu", 29},
	{"Db", 105},
	{"Ds", 110},
	{"Dy", 66},
	{"Er", 68},
	{"Es", 99},
	{"Eu", 63},
	{"F", 9},
	{"Fe", 26},
	{"Fm", 100},
	{"Fr", 87},
	{"Ga", 31},
	{"Gd", 64},
	{"Ge", 32},
	{"H", 1},
	{"He", 2},
	{"Hf", 72},
	{"Hg", 80},
	{"Ho", 67},
	{"Hs", 108},
	{"I", 53},
	{"In", 49},
	{"Ir", 77},
	{"K", 19},
	{"Kr", 36},
	{"La", 57},
	{"Li", 3},
	{"Lr", 103},
	{"Lu", 71},
	{"Md", 101},
	{"Mg", 12},
	{"Mn", 25},
	{"Mo", 42},
	{"Mt", 109},
	{"N", 7},
	{"Na", 11},
	{"Nb", 41},
	{"Nd", 60},
	{"Ne", 10},
	{"Ni", 28},
	{"No", 102},
	{"Np", 93},
	{"O", 8},
	{"Os", 76},
	{"P", 15},
	{"Pa", 91},
	{"Pb", 82},
	{"Pd", 46},
	{"Pm", 61},
	{"Po", 84},
	{"Pr", 59},
	{"Pt", 78},
	{"Pu", 94},
	{"Ra", 88},
	{"Rb", 37},
	{"Re", 75},
	{"Rf", 104},
	{"Rg", 111},
	{"Rh", 45},
	{"Rn", 86},
	{"Ru", 44},
	{"S", 16},
	{"Sb", 51},
	{"Sc", 21},
	{"Se", 34},
	{"Sg", 106},
	{"Si", 14},
	{"Sm", 62},
	{"Sn", 50},
	{"Sr", 38},
	{"Ta", 73},
	{"Tb", 65},
	{"Tc", 43},
	{"Te", 52},
	{"Th", 90},
	{"Ti", 22},
	{"Tl", 81},
	{"Tm", 69},
	{"U", 92},
	{"Uuh", 116},
	{"Uuo", 118},
	{"Uup", 115},
	{"Uuq", 114},
	{"Uus", 117},
	{"Uut", 113},
	{"V", 23},
	{"W", 74},
	{"Xe", 54},
	{"Y", 39},
	{"Yb", 70},
	{"Zn", 30},
	{"Zr", 40}
};

map<string,int> atomicSymbolMap(atomicSymbolArray, atomicSymbolArray+nElements);
