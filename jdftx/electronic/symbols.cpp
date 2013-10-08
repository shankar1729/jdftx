/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <electronic/symbols.h>

EnumStringMap<AtomicSymbol> atomicSymbolMap
(	AtomicSymbol::H, "H",
	AtomicSymbol::He, "He",
	AtomicSymbol::Li, "Li",
	AtomicSymbol::Be, "Be",
	AtomicSymbol::B, "B",
	AtomicSymbol::C, "C",
	AtomicSymbol::N, "N",
	AtomicSymbol::O, "O",
	AtomicSymbol::F, "F",
	AtomicSymbol::Ne, "Ne",
	AtomicSymbol::Na, "Na",
	AtomicSymbol::Mg, "Mg",
	AtomicSymbol::Al, "Al",
	AtomicSymbol::Si, "Si",
	AtomicSymbol::P, "P",
	AtomicSymbol::S, "S",
	AtomicSymbol::Cl, "Cl",
	AtomicSymbol::Ar, "Ar",
	AtomicSymbol::K, "K",
	AtomicSymbol::Ca, "Ca",
	AtomicSymbol::Sc, "Sc",
	AtomicSymbol::Ti, "Ti",
	AtomicSymbol::V, "V",
	AtomicSymbol::Cr, "Cr",
	AtomicSymbol::Mn, "Mn",
	AtomicSymbol::Fe, "Fe",
	AtomicSymbol::Co, "Co",
	AtomicSymbol::Ni, "Ni",
	AtomicSymbol::Cu, "Cu",
	AtomicSymbol::Zn, "Zn",
	AtomicSymbol::Ga, "Ga",
	AtomicSymbol::Ge, "Ge",
	AtomicSymbol::As, "As",
	AtomicSymbol::Se, "Se",
	AtomicSymbol::Br, "Br",
	AtomicSymbol::Kr, "Kr",
	AtomicSymbol::Rb, "Rb",
	AtomicSymbol::Sr, "Sr",
	AtomicSymbol::Y, "Y",
	AtomicSymbol::Zr, "Zr",
	AtomicSymbol::Nb, "Nb",
	AtomicSymbol::Mo, "Mo",
	AtomicSymbol::Tc, "Tc",
	AtomicSymbol::Ru, "Ru",
	AtomicSymbol::Rh, "Rh",
	AtomicSymbol::Pd, "Pd",
	AtomicSymbol::Ag, "Ag",
	AtomicSymbol::Cd, "Cd",
	AtomicSymbol::In, "In",
	AtomicSymbol::Sn, "Sn",
	AtomicSymbol::Sb, "Sb",
	AtomicSymbol::Te, "Te",
	AtomicSymbol::I, "I",
	AtomicSymbol::Xe, "Xe",
	AtomicSymbol::Cs, "Cs",
	AtomicSymbol::Ba, "Ba",
	AtomicSymbol::La, "La",
	AtomicSymbol::Ce, "Ce",
	AtomicSymbol::Pr, "Pr",
	AtomicSymbol::Nd, "Nd",
	AtomicSymbol::Pm, "Pm",
	AtomicSymbol::Sm, "Sm",
	AtomicSymbol::Eu, "Eu",
	AtomicSymbol::Gd, "Gd",
	AtomicSymbol::Tb, "Tb",
	AtomicSymbol::Dy, "Dy",
	AtomicSymbol::Ho, "Ho",
	AtomicSymbol::Er, "Er",
	AtomicSymbol::Tm, "Tm",
	AtomicSymbol::Yb, "Yb",
	AtomicSymbol::Lu, "Lu",
	AtomicSymbol::Hf, "Hf",
	AtomicSymbol::Ta, "Ta",
	AtomicSymbol::W, "W",
	AtomicSymbol::Re, "Re",
	AtomicSymbol::Os, "Os",
	AtomicSymbol::Ir, "Ir",
	AtomicSymbol::Pt, "Pt",
	AtomicSymbol::Au, "Au",
	AtomicSymbol::Hg, "Hg",
	AtomicSymbol::Tl, "Tl",
	AtomicSymbol::Pb, "Pb",
	AtomicSymbol::Bi, "Bi",
	AtomicSymbol::Po, "Po",
	AtomicSymbol::At, "At",
	AtomicSymbol::Rn, "Rn",
	AtomicSymbol::Fr, "Fr",
	AtomicSymbol::Ra, "Ra",
	AtomicSymbol::Ac, "Ac",
	AtomicSymbol::Th, "Th",
	AtomicSymbol::Pa, "Pa",
	AtomicSymbol::U, "U",
	AtomicSymbol::Np, "Np",
	AtomicSymbol::Pu, "Pu",
	AtomicSymbol::Am, "Am",
	AtomicSymbol::Cm, "Cm",
	AtomicSymbol::Bk, "Bk",
	AtomicSymbol::Cf, "Cf",
	AtomicSymbol::Es, "Es",
	AtomicSymbol::Fm, "Fm",
	AtomicSymbol::Md, "Md",
	AtomicSymbol::No, "No",
	AtomicSymbol::Lr, "Lr",
	AtomicSymbol::Rf, "Rf",
	AtomicSymbol::Db, "Db",
	AtomicSymbol::Sg, "Sg",
	AtomicSymbol::Bh, "Bh",
	AtomicSymbol::Hs, "Hs",
	AtomicSymbol::Mt, "Mt",
	AtomicSymbol::Ds, "Ds",
	AtomicSymbol::Rg, "Rg",
	AtomicSymbol::Cn, "Cn",
	AtomicSymbol::Uut, "Uut",
	AtomicSymbol::Uuq, "Uuq",
	AtomicSymbol::Uup, "Uup",
	AtomicSymbol::Uuh, "Uuh",
	AtomicSymbol::Uus, "Uus",
	AtomicSymbol::Uuo, "Uuo"
);


double atomicMass(AtomicSymbol symbol)
{	if(symbol==AtomicSymbol::Ac) return 227.028000;
	if(symbol==AtomicSymbol::Al) return 26.981539;
	if(symbol==AtomicSymbol::Am) return 243.000000;
	if(symbol==AtomicSymbol::Sb) return 121.757000;
	if(symbol==AtomicSymbol::Ar) return 39.948000;
	if(symbol==AtomicSymbol::As) return 74.921590;
	if(symbol==AtomicSymbol::At) return 210.000000;
	if(symbol==AtomicSymbol::Ba) return 137.327000;
	if(symbol==AtomicSymbol::Bk) return 247.000000;
	if(symbol==AtomicSymbol::Be) return 9.012182;
	if(symbol==AtomicSymbol::Bi) return 208.980370;
	if(symbol==AtomicSymbol::Bh) return 262.000000;
	if(symbol==AtomicSymbol::B) return 10.811000;
	if(symbol==AtomicSymbol::Br) return 79.904000;
	if(symbol==AtomicSymbol::Cd) return 112.411000;
	if(symbol==AtomicSymbol::Ca) return 40.078000;
	if(symbol==AtomicSymbol::Cf) return 251.000000;
	if(symbol==AtomicSymbol::C) return 12.011000;
	if(symbol==AtomicSymbol::Ce) return 140.115000;
	if(symbol==AtomicSymbol::Cs) return 132.905430;
	if(symbol==AtomicSymbol::Cl) return 35.452700;
	if(symbol==AtomicSymbol::Cr) return 51.996100;
	if(symbol==AtomicSymbol::Co) return 58.933200;
	if(symbol==AtomicSymbol::Cu) return 63.546000;
	if(symbol==AtomicSymbol::Cm) return 247.000000;
	if(symbol==AtomicSymbol::Db) return 262.000000;
	if(symbol==AtomicSymbol::Dy) return 162.500000;
	if(symbol==AtomicSymbol::Es) return 252.000000;
	if(symbol==AtomicSymbol::Er) return 167.260000;
	if(symbol==AtomicSymbol::Eu) return 151.965000;
	if(symbol==AtomicSymbol::Fm) return 257.000000;
	if(symbol==AtomicSymbol::F) return 18.998403;
	if(symbol==AtomicSymbol::Fr) return 223.000000;
	if(symbol==AtomicSymbol::Gd) return 157.250000;
	if(symbol==AtomicSymbol::Ga) return 69.723000;
	if(symbol==AtomicSymbol::Ge) return 72.610000;
	if(symbol==AtomicSymbol::Au) return 196.966540;
	if(symbol==AtomicSymbol::Hf) return 178.490000;
	if(symbol==AtomicSymbol::Hs) return 265.000000;
	if(symbol==AtomicSymbol::He) return 4.002602;
	if(symbol==AtomicSymbol::Ho) return 164.930320;
	if(symbol==AtomicSymbol::H) return 1.007940;
	if(symbol==AtomicSymbol::In) return 114.820000;
	if(symbol==AtomicSymbol::I) return 126.904470;
	if(symbol==AtomicSymbol::Ir) return 192.220000;
	if(symbol==AtomicSymbol::Fe) return 55.847000;
	if(symbol==AtomicSymbol::Kr) return 83.800000;
	if(symbol==AtomicSymbol::La) return 138.905500;
	if(symbol==AtomicSymbol::Lr) return 262.000000;
	if(symbol==AtomicSymbol::Pb) return 207.200000;
	if(symbol==AtomicSymbol::Li) return 6.941000;
	if(symbol==AtomicSymbol::Lu) return 174.967000;
	if(symbol==AtomicSymbol::Mg) return 24.305000;
	if(symbol==AtomicSymbol::Mn) return 54.938050;
	if(symbol==AtomicSymbol::Mt) return 266.000000;
	if(symbol==AtomicSymbol::Md) return 258.000000;
	if(symbol==AtomicSymbol::Hg) return 200.590000;
	if(symbol==AtomicSymbol::Mo) return 95.940000;
	if(symbol==AtomicSymbol::Nd) return 144.240000;
	if(symbol==AtomicSymbol::Ne) return 20.179700;
	if(symbol==AtomicSymbol::Np) return 237.048000;
	if(symbol==AtomicSymbol::Ni) return 58.693400;
	if(symbol==AtomicSymbol::Nb) return 92.906380;
	if(symbol==AtomicSymbol::N) return 14.006740;
	if(symbol==AtomicSymbol::No) return 259.000000;
	if(symbol==AtomicSymbol::Os) return 190.200000;
	if(symbol==AtomicSymbol::O) return 15.999400;
	if(symbol==AtomicSymbol::Pd) return 106.420000;
	if(symbol==AtomicSymbol::P) return 30.973762;
	if(symbol==AtomicSymbol::Pt) return 195.080000;
	if(symbol==AtomicSymbol::Pu) return 244.000000;
	if(symbol==AtomicSymbol::Po) return 209.000000;
	if(symbol==AtomicSymbol::K) return 39.098300;
	if(symbol==AtomicSymbol::Pr) return 140.907650;
	if(symbol==AtomicSymbol::Pm) return 145.000000;
	if(symbol==AtomicSymbol::Pa) return 231.035900;
	if(symbol==AtomicSymbol::Ra) return 226.025000;
	if(symbol==AtomicSymbol::Rn) return 222.000000;
	if(symbol==AtomicSymbol::Re) return 186.207000;
	if(symbol==AtomicSymbol::Rh) return 102.905500;
	if(symbol==AtomicSymbol::Rb) return 85.467800;
	if(symbol==AtomicSymbol::Ru) return 101.070000;
	if(symbol==AtomicSymbol::Rf) return 261.000000;
	if(symbol==AtomicSymbol::Sm) return 150.360000;
	if(symbol==AtomicSymbol::Sc) return 44.955910;
	if(symbol==AtomicSymbol::Sg) return 263.000000;
	if(symbol==AtomicSymbol::Se) return 78.960000;
	if(symbol==AtomicSymbol::Si) return 28.085500;
	if(symbol==AtomicSymbol::Ag) return 107.868200;
	if(symbol==AtomicSymbol::Na) return 22.989768;
	if(symbol==AtomicSymbol::Sr) return 87.620000;
	if(symbol==AtomicSymbol::S) return 32.066000;
	if(symbol==AtomicSymbol::Ta) return 180.947900;
	if(symbol==AtomicSymbol::Tc) return 98.000000;
	if(symbol==AtomicSymbol::Te) return 127.600000;
	if(symbol==AtomicSymbol::Tb) return 158.925340;
	if(symbol==AtomicSymbol::Tl) return 204.383300;
	if(symbol==AtomicSymbol::Th) return 232.038100;
	if(symbol==AtomicSymbol::Tm) return 168.934210;
	if(symbol==AtomicSymbol::Sn) return 118.710000;
	if(symbol==AtomicSymbol::Ti) return 47.880000;
	if(symbol==AtomicSymbol::W) return 183.850000;
	if(symbol==AtomicSymbol::U) return 238.028900;
	if(symbol==AtomicSymbol::V) return 50.941500;
	if(symbol==AtomicSymbol::Xe) return 131.290000;
	if(symbol==AtomicSymbol::Yb) return 173.040000;
	if(symbol==AtomicSymbol::Y) return 88.905850;
	if(symbol==AtomicSymbol::Zn) return 65.390000;
	if(symbol==AtomicSymbol::Zr) return 91.224000;
	die("Atomic mass unavilable for element %s\n", atomicSymbolMap.getString(symbol));
}
