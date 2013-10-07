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
{	switch(symbol)
	{	case AtomicSymbol::Ac: return 227.028000;
		case AtomicSymbol::Al: return 26.981539;
		case AtomicSymbol::Am: return 243.000000;
		case AtomicSymbol::Sb: return 121.757000;
		case AtomicSymbol::Ar: return 39.948000;
		case AtomicSymbol::As: return 74.921590;
		case AtomicSymbol::At: return 210.000000;
		case AtomicSymbol::Ba: return 137.327000;
		case AtomicSymbol::Bk: return 247.000000;
		case AtomicSymbol::Be: return 9.012182;
		case AtomicSymbol::Bi: return 208.980370;
		case AtomicSymbol::Bh: return 262.000000;
		case AtomicSymbol::B: return 10.811000;
		case AtomicSymbol::Br: return 79.904000;
		case AtomicSymbol::Cd: return 112.411000;
		case AtomicSymbol::Ca: return 40.078000;
		case AtomicSymbol::Cf: return 251.000000;
		case AtomicSymbol::C: return 12.011000;
		case AtomicSymbol::Ce: return 140.115000;
		case AtomicSymbol::Cs: return 132.905430;
		case AtomicSymbol::Cl: return 35.452700;
		case AtomicSymbol::Cr: return 51.996100;
		case AtomicSymbol::Co: return 58.933200;
		case AtomicSymbol::Cu: return 63.546000;
		case AtomicSymbol::Cm: return 247.000000;
		case AtomicSymbol::Db: return 262.000000;
		case AtomicSymbol::Dy: return 162.500000;
		case AtomicSymbol::Es: return 252.000000;
		case AtomicSymbol::Er: return 167.260000;
		case AtomicSymbol::Eu: return 151.965000;
		case AtomicSymbol::Fm: return 257.000000;
		case AtomicSymbol::F: return 18.998403;
		case AtomicSymbol::Fr: return 223.000000;
		case AtomicSymbol::Gd: return 157.250000;
		case AtomicSymbol::Ga: return 69.723000;
		case AtomicSymbol::Ge: return 72.610000;
		case AtomicSymbol::Au: return 196.966540;
		case AtomicSymbol::Hf: return 178.490000;
		case AtomicSymbol::Hs: return 265.000000;
		case AtomicSymbol::He: return 4.002602;
		case AtomicSymbol::Ho: return 164.930320;
		case AtomicSymbol::H: return 1.007940;
		case AtomicSymbol::In: return 114.820000;
		case AtomicSymbol::I: return 126.904470;
		case AtomicSymbol::Ir: return 192.220000;
		case AtomicSymbol::Fe: return 55.847000;
		case AtomicSymbol::Kr: return 83.800000;
		case AtomicSymbol::La: return 138.905500;
		case AtomicSymbol::Lr: return 262.000000;
		case AtomicSymbol::Pb: return 207.200000;
		case AtomicSymbol::Li: return 6.941000;
		case AtomicSymbol::Lu: return 174.967000;
		case AtomicSymbol::Mg: return 24.305000;
		case AtomicSymbol::Mn: return 54.938050;
		case AtomicSymbol::Mt: return 266.000000;
		case AtomicSymbol::Md: return 258.000000;
		case AtomicSymbol::Hg: return 200.590000;
		case AtomicSymbol::Mo: return 95.940000;
		case AtomicSymbol::Nd: return 144.240000;
		case AtomicSymbol::Ne: return 20.179700;
		case AtomicSymbol::Np: return 237.048000;
		case AtomicSymbol::Ni: return 58.693400;
		case AtomicSymbol::Nb: return 92.906380;
		case AtomicSymbol::N: return 14.006740;
		case AtomicSymbol::No: return 259.000000;
		case AtomicSymbol::Os: return 190.200000;
		case AtomicSymbol::O: return 15.999400;
		case AtomicSymbol::Pd: return 106.420000;
		case AtomicSymbol::P: return 30.973762;
		case AtomicSymbol::Pt: return 195.080000;
		case AtomicSymbol::Pu: return 244.000000;
		case AtomicSymbol::Po: return 209.000000;
		case AtomicSymbol::K: return 39.098300;
		case AtomicSymbol::Pr: return 140.907650;
		case AtomicSymbol::Pm: return 145.000000;
		case AtomicSymbol::Pa: return 231.035900;
		case AtomicSymbol::Ra: return 226.025000;
		case AtomicSymbol::Rn: return 222.000000;
		case AtomicSymbol::Re: return 186.207000;
		case AtomicSymbol::Rh: return 102.905500;
		case AtomicSymbol::Rb: return 85.467800;
		case AtomicSymbol::Ru: return 101.070000;
		case AtomicSymbol::Rf: return 261.000000;
		case AtomicSymbol::Sm: return 150.360000;
		case AtomicSymbol::Sc: return 44.955910;
		case AtomicSymbol::Sg: return 263.000000;
		case AtomicSymbol::Se: return 78.960000;
		case AtomicSymbol::Si: return 28.085500;
		case AtomicSymbol::Ag: return 107.868200;
		case AtomicSymbol::Na: return 22.989768;
		case AtomicSymbol::Sr: return 87.620000;
		case AtomicSymbol::S: return 32.066000;
		case AtomicSymbol::Ta: return 180.947900;
		case AtomicSymbol::Tc: return 98.000000;
		case AtomicSymbol::Te: return 127.600000;
		case AtomicSymbol::Tb: return 158.925340;
		case AtomicSymbol::Tl: return 204.383300;
		case AtomicSymbol::Th: return 232.038100;
		case AtomicSymbol::Tm: return 168.934210;
		case AtomicSymbol::Sn: return 118.710000;
		case AtomicSymbol::Ti: return 47.880000;
		case AtomicSymbol::W: return 183.850000;
		case AtomicSymbol::U: return 238.028900;
		case AtomicSymbol::V: return 50.941500;
		case AtomicSymbol::Xe: return 131.290000;
		case AtomicSymbol::Yb: return 173.040000;
		case AtomicSymbol::Y: return 88.905850;
		case AtomicSymbol::Zn: return 65.390000;
		case AtomicSymbol::Zr: return 91.224000;
		default: die("Atomic mass unavilable for element %s\n", atomicSymbolMap.getString(symbol));
	}
}
