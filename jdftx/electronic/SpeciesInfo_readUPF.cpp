/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <electronic/SpeciesInfo.h>
#include <electronic/Everything.h>
#include <electronic/matrix.h>
#include <electronic/symbols.h>
#include "SpeciesInfo_internal.h"

//! Read an XML tag from the specified stream and parse any key-value pairs 
struct XMLtag
{	istream& is;
	string name;
	std::map<string,string> attributes; //!< list of key-value pairs
	bool closed; //!< whether the tag has already been closed (either by the empty /> termination or by the matching closing tag)
	
	XMLtag(istream& is) : is(is), closed(false)
	{	//Read tag name:
		while(true)
		{	char cStart = getNonSpace();
			if(cStart != '<')
				die("  XML parse error: expecting tag starting with '<'; found '%c%s' instead.\n", cStart, readToken().c_str());
			if(is.peek()=='/') //Found a closing tag instead of an opening one
			{	is.putback('<');
				name.clear(); //indicate that the tag is null
				closed = true;
				return;
			}
			name = readToken();
			//Check for and ignore comments:
			if(name == "!--")
			{	char c1 = ' ', c2 = ' ';
				while(true)
				{	char c = is.get();
					if(c2 == '-' && c1 == '-' && c == '>')
						break;
					if(is.eof()) die("  XML parse error: unterminated comment.\n");
					c2 = c1; c1 = c;
				}
			}
			else break;
		}
		//Read attributes:
		while(true)
		{	//Check for closing:
			char c = getNonSpace();
			if(c=='>') break; //End of attributes
			else if(c=='/')
			{	if(is.get()!='>') die("  XML parse error in tag '%s': '/' not followed by '>'\n", name.c_str());
				closed = true; //empty tag (no data section)
				break;
			}
			else is.putback(c);
			//Get attribute name:
			string key = readToken();
			//Get separator:
			c = getNonSpace();
			if(c != '=') die("  XML parse error in tag '%s': first non-space character after attribute name must be '='; found '%c' instead\n", name.c_str(), c);
			//Get value:
			char cOpen = getNonSpace();
			if(!(cOpen=='"' || cOpen=='\''))
				die("  XML parse error in tag '%s': attribute value must start with a single or double quote; found '%c' instead\n", name.c_str(), cOpen);
			string value;
			while((c=is.get()) != cOpen)
			{	value.push_back(c); //read upto terminating quote
				if(is.eof()) die("  XML parse error in tag '%s': file ended with unterminated string\n", name.c_str());
			}
			//Save key value pair:
			attributes[key] = value;
		}
	}
	
	void close() //!< if not already closed, read the closing tag
	{	if(!closed)
		{	char c = getNonSpace(); if(c != '<') die("  XML parse error: expecting '</%s'; found '%c' instead\n", name.c_str(), c);
			c = is.get();  if(c != '/') die("  XML parse error: expecting '</%s'; found '<%c' instead\n", name.c_str(), c);
			string tag = readToken();
			if(tag != name) die("  XML parse error: expecting '</%s'; found '</%s' instead\n", name.c_str(), tag.c_str());
			c = getNonSpace(); if(c != '>') die("  XML parse error: expecting '>' to close '</%s'; found '%c' instead\n", name.c_str(), c);
			closed = true;
		}
	}
	
	void ignoreAndClose() //!< Ignore any data or sub-tags till matching close tag is found
	{	if(!closed)
		{	while(true)
			{	char c = is.get();
				if(c == '<')
				{	c = is.get();
					if(c == '/')
					{	string tag = readToken();
						if(tag == name)
						{	c = getNonSpace(); if(c != '>') die("  XML parse error: expecting '>' to close '</%s'; found '%c' instead\n", name.c_str(), c);
							closed = true;
							return;
						}
					}
				}
				if(is.eof()) die("  XML parse error: file ended prematurely while searching for closing tag '</%s'\n", name.c_str());
			}
		}
	}
	
	string getAttribute(string key, bool required=true) //!< Wrapper to access attributes with error message
	{	auto iter = attributes.find(key);
		if(iter == attributes.end())
		{	if(required) die("  Compulsory attribute '%s' missing in section '%s'\n", key.c_str(), name.c_str());
			return string();
		}
		return iter->second;
	}
	
	std::vector<double> readData(size_t nElem)
	{	std::vector<double> out(nElem);
		for(double& x: out)
		{	is >> x;
			if(is.eof()) die("  XML parse error: file ended prematurely while reading data in section '%s'\n", name.c_str());
		}
		return out;
	}
	
private:
	//! Get first character which is not whitespace
	char getNonSpace()
	{	char c=is.get();
		while(isspace(c)) c=is.get();
		return c;
	}
	
	//! Read a valid XML tag or attribute name
	string readToken()
	{	string token;
		char c = getNonSpace();
		while(true)
		{	if(isspace(c) || c=='=' || c=='>' || c=='/' || c=='<')
			{	is.putback(c);
				break;
			}
			token.push_back(c);
			c = is.get();
		}
		return token;
	}
};

void SpeciesInfo::readUPF(istream& is)
{	
	int lMax = 0; //max angular momentum
	int nGrid = 0; //number of grid points in radial mesh
	int nBeta = 0; //number of projectors
	int nPsi = 0; //number of atomic orbitals
	std::vector<double> rGrid, drGrid; //radial grid and integration factor
	std::set<double> coreRadii; //ordered list of various cutoff radii
	std::vector<int> lNL; std::vector<double> jNL; //orbital and total angular momentum per projector
	std::vector<int> lPsi; std::vector<double> jPsi; //orbital and total angular momentum per orbital
	
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	int nGridNL = int(ceil(e->gInfo.GmaxSphere/dG))+5;

	//Read XML file:
	XMLtag tagUPF(is);
	if(tagUPF.name != "UPF")
		die("  First XML tag must be 'UPF'; found '%s' instead.\n"
			"  Perhaps this is an old-format UPF file; if so, use upf2upf2\n"
			"  distributed with Quantum Espresso to convert it to UPF >=2.0.1.\n", tagUPF.name.c_str());
	while(true)
	{	XMLtag tag(is);
		if(!tag.name.length()) break;
		if(tag.name == "PP_INFO" || tag.name == "PP_RHOATOM")
		{	tag.ignoreAndClose(); //Ignore silently
		}
		else if(tag.name == "PP_HEADER")
		{	//Species identity:
			string element = tag.getAttribute("element"); trim(element);
			AtomicSymbol atSym = AtomicSymbol::H;
			if(!atomicSymbolMap.getEnum(element.c_str(), atSym))
				die("  Could not determine atomic number for element '%s'.\n"
					"  Please edit pseudopotential to use the standard chemical symbol.\n", element.c_str());
			atomicNumber = int(atSym);
			logPrintf("  '%s' pseudopotential, '%s' functional\n", element.c_str(), tag.getAttribute("functional").c_str());
			//Non essential info:
			string generated = tag.getAttribute("generated", false); if(generated.length()) logPrintf("  %s\n", generated.c_str());
			string comment = tag.getAttribute("comment", false); if(comment.length()) logPrintf("  %s\n", comment.c_str());
			string author = tag.getAttribute("author", false); if(author.length()) logPrintf("  Author: %s", author.c_str());
			string date = tag.getAttribute("date", false); if(date.length()) logPrintf("  Date: %s", date.c_str());
			if(author.length() || date.length()) logPrintf(".\n");
			//Check for unsupported types:
			if(tag.getAttribute("is_paw") == "T")
				die("  PAW datasets are not yet supported.\n");
			//Valence properties:
			Z = atof(tag.getAttribute("z_valence").c_str());
			lMax = atoi(tag.getAttribute("l_max").c_str());
			nGrid = atoi(tag.getAttribute("mesh_size").c_str());
			nBeta = atoi(tag.getAttribute("number_of_proj").c_str());
			nPsi = atoi(tag.getAttribute("number_of_wfc").c_str());
			logPrintf("  %lg valence electrons, %d orbitals, %d projectors, %d radial grid points, with lMax = %d\n", Z, nPsi, nBeta, nGrid, lMax);
		}
		else if(tag.name == "PP_MESH")
		{	while(true)
			{	XMLtag tagMesh(is);
				if(!tagMesh.name.length()) break;
				if(tagMesh.name == "PP_R") rGrid = tagMesh.readData(nGrid);
				else if(tagMesh.name == "PP_RAB") drGrid = tagMesh.readData(nGrid);
				else
				{	logPrintf("  NOTE: ignored section '%s'\n", tagMesh.name.c_str());
					tagMesh.ignoreAndClose();
				}
				tagMesh.close();
			}
		}
		else if(tag.name == "PP_NLCC")
		{	RadialFunctionR nCore(nGrid);
			nCore.set(rGrid, drGrid);
			nCore.f = tag.readData(nGrid);
			setCore(nCore);
		}
		else if(tag.name == "PP_LOCAL")
		{	RadialFunctionR Vloc(nGrid);
			Vloc.set(rGrid, drGrid);
			Vloc.f = tag.readData(nGrid);
			for(int i=0; i<nGrid; i++)
				Vloc.f[i] = 0.5*Vloc.f[i] + Z*(rGrid[i] ? 1./rGrid[i] : 0); //Convert from Ry to Eh and remove Z/r part
			logPrintf("  Transforming local potential to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridLoc);
			Vloc.transform(0, dG, nGridLoc, VlocRadial);
		}
		else if(tag.name == "PP_NONLOCAL")
		{	lNL.assign(nBeta, -1); //angular momentum per projector
			std::vector<RadialFunctionR> Vnl(nBeta, RadialFunctionR(nGrid)); //projectors
			std::vector<std::vector<double> > D(nBeta); //descreened D matrix (Mnl in our notation)
			std::vector<std::vector<double> > Q(nBeta); //overlap augmentation matrix (goes to Qint in SpeciesInfo)
			std::vector<std::vector<RadialFunctionR> > Qr(nBeta, std::vector<RadialFunctionR>(nBeta, RadialFunctionR(nGrid))); //overlap radial function
			std::vector<std::vector<std::vector<RadialFunctionR> > > Qlr(nBeta, std::vector<std::vector<RadialFunctionR> >(nBeta)); //explicitly l-dependent overlap radial function
			std::vector<std::vector<std::vector<std::vector<double> > > > Qcoeff(nBeta, std::vector<std::vector<std::vector<double> > >(nBeta));
			std::vector<double> rInner; //Replace Qr / Qlr with Taylor expansion from Qcoeff within rInner
			bool q_with_l = false; //File contains Qlr if true and Qr if false (see above)
			
			while(true)
			{	XMLtag tagNL(is);
				if(!tagNL.name.length()) break;
				if(tagNL.name.substr(0,8) == "PP_BETA.") //K-B projector
				{	int iBeta = atoi(tagNL.name.substr(8).c_str());
					if(iBeta<1 || iBeta>nBeta)
						die("  Invalid projector index %d (not in range [1,%d]) in section '%s'\n", iBeta, nBeta, tagNL.name.c_str());
					iBeta--; //to zero-based index
					lNL[iBeta] = atoi(tagNL.getAttribute("angular_momentum").c_str());
					if(lNL[iBeta]<0 || lNL[iBeta]>lMax)
						die("  Invalid projector angular momentum %d (not in range [0, lMax=%d]) in section '%s'\n", lNL[iBeta], lMax, tagNL.name.c_str());
					Vnl[iBeta].set(rGrid, drGrid);
					Vnl[iBeta].f = tagNL.readData(nGrid);
				}
				else if(tagNL.name == "PP_DIJ")
				{	for(std::vector<double>& row: D)
						row =tagNL.readData(nBeta);
				}
				else if(tagNL.name == "PP_AUGMENTATION")
				{	q_with_l = (tagNL.getAttribute("q_with_l")=="T");
					int nCoeff = atoi(tagNL.getAttribute("nqf").c_str());
					int nQLC = atoi(tagNL.getAttribute("nqlc").c_str()); if(!nQLC) nQLC = 2*lMax+1;
					while(true)
					{	XMLtag tagAug(is);
						if(!tagAug.name.length()) break;
						if(tagAug.name == "PP_Q")
						{	for(std::vector<double>& row: Q)
								row = tagAug.readData(nBeta);
						}
						else if(tagAug.name == "PP_QFCOEF")
						{	for(auto& Qcoeff_i: Qcoeff)
								for(auto& Qcoeff_ij: Qcoeff_i)
								{	Qcoeff_ij.resize(nQLC);
									for(auto& Qcoeff_ijl: Qcoeff_ij)
										Qcoeff_ijl =  tagAug.readData(nCoeff);
								}
						}
						else if(tagAug.name == "PP_RINNER")
						{	rInner = tagAug.readData(nQLC);
						}
						else if(tagAug.name.substr(0,6)=="PP_QIJ")
						{	istringstream iss(tagAug.name);
							string buf; getline(iss, buf, '.');
							getline(iss, buf, '.'); int i = atoi(buf.c_str());
							getline(iss, buf, '.'); int j = atoi(buf.c_str());
							getline(iss, buf, '.'); int l = atoi(buf.c_str()); //this will cause a harmless eof on iss when q_with_l = false
							if(i<1 || i>nBeta) die("  Invalid projector index %d (not in range [1,%d]) in section '%s'\n", i, nBeta, tagAug.name.c_str());
							if(j<1 || j>nBeta) die("  Invalid projector index %d (not in range [1,%d]) in section '%s'\n", j, nBeta, tagAug.name.c_str());
							if(q_with_l && (l<0 || l>2*lMax))
								die("  Invalid augmentation angular momentum %d (not in range [0, 2*lMax=%d]) in section '%s'\n", l, 2*lMax, tagNL.name.c_str());
							i--; j--; //to zero-based indices
							if(i>j) std::swap(i,j); //store results only for i<=j (symmetry)
							RadialFunctionR* Qtarget = 0;
							if(q_with_l)
							{	if(!Qlr[i][j].size()) Qlr[i][j].resize(2*lMax+1, RadialFunctionR(nGrid));
								Qtarget = &Qlr[i][j][l];
							}
							else Qtarget = &Qr[i][j];
							Qtarget->set(rGrid, drGrid);
							Qtarget->f = tagAug.readData(nGrid);
						}
						else
						{	logPrintf("  NOTE: ignored section '%s'\n", tagAug.name.c_str());
							tagAug.ignoreAndClose();
						}
						tagAug.close();
					}
				}
				else
				{	logPrintf("  NOTE: ignored section '%s'\n", tagNL.name.c_str());
					tagNL.ignoreAndClose();
				}
				tagNL.close();
			}
			
			if(nBeta > 0)
			{	if(!D[0].size()) die("  Nonlocal pseudopotential 'D' matrix has not defined.\n");
				//Check and transform projectors:
				logPrintf("  Transforming nonlocal projectors to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridNL);
				VnlRadial.resize(lMax+1);
				std::vector<std::vector<int> > lBeta(lMax+1); //index map from l and p to combined index iBeta
				for(int iBeta=0; iBeta<nBeta; iBeta++)
				{	int l = lNL[iBeta];
					if(l<0) die("  Nonlocal projector number %d has not been defined.\n", iBeta+1);
					lBeta[l].push_back(iBeta);
					VnlRadial[l].push_back(RadialFunctionG());
					Vnl[iBeta].set(rGrid, drGrid);
					for(int i=0; i<nGrid; i++)
						Vnl[iBeta].f[i] *= (rGrid[i] ? 1./rGrid[i] : 0);
					Vnl[iBeta].transform(l, dG, nGridNL, VnlRadial[l].back());
					//Determine core radius:
					for(int i=nGrid; i>=0; i--)
						if(4*M_PI*rGrid[i]*rGrid[i]*drGrid[i] * fabs(D[iBeta][iBeta]) * Vnl[iBeta].f[i]*Vnl[iBeta].f[i] > 1e-3)
						{	coreRadii.insert(rGrid[i]);
							break;
						}
				}
				//Set Mnl:
				Mnl.resize(lMax+1);
				for(int l=0; l<=lMax; l++)
				{	int nProj = lBeta[l].size();
					Mnl[l].init(nProj, nProj);
					for(int p1=0; p1<nProj; p1++)
						for(int p2=0; p2<nProj; p2++)
						{	int iBeta = lBeta[l][p1];
							int jBeta = lBeta[l][p2];
							if(iBeta>jBeta) std::swap(iBeta,jBeta);
							Mnl[l].data()[Mnl[l].index(p1,p2)] = D[iBeta][jBeta]*0.5; //Convert to Eh
						}
				}
				//Q radial functions and integrals:
				if(Q[0].size())
				{	logPrintf("  Transforming density augmentations to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridLoc);
					Qint.resize(lMax+1);
					for(int l=0; l<=lMax; l++) if(lBeta[l].size())
						Qint[l] = zeroes(lBeta[l].size(), lBeta[l].size());
					for(int l1=0; l1<=lMax; l1++) for(int p1=0; p1<int(lBeta[l1].size()); p1++)
					{	int iBeta = lBeta[l1][p1];
						for(int l2=0; l2<=lMax; l2++) for(int p2=0; p2<int(lBeta[l2].size()); p2++)
						{	int jBeta = lBeta[l2][p2];
							if(iBeta<=jBeta) //store and use upper triangular only
							{	for(int l=abs(l1-l2); l<=l1+l2; l+=2)
								{	const RadialFunctionR& Qsrc = q_with_l ? Qlr[iBeta][jBeta][l] : Qr[iBeta][jBeta];
									const std::vector<double>* coeff = (rInner.size() && Qcoeff[iBeta][jBeta].size()) ?  &Qcoeff[iBeta][jBeta][l] : 0;
									//Replace the inner section with the l-dependent pseudization:
									RadialFunctionR Qijl(nGrid); Qijl.set(rGrid, drGrid);
									bool isNonzero = false; int iLastNZ=0; double Qabs = 0.;
									for(int i=0; i<nGrid; i++)
									{	if(coeff && rGrid[i]<rInner[l])
										{	double val=0.0, rSq=rGrid[i]*rGrid[i], rPow=pow(rGrid[i],l);
											for(auto c: *coeff)
											{	val += c*rPow;
												rPow *= rSq;
											}
											Qijl.f[i] = val;
										}
										else Qijl.f[i] = Qsrc.f[i]/pow(rGrid[i],2);
										//Radial extent check:
										double dQabs = fabs(Qijl.f[i]) * rGrid[i]*rGrid[i] * drGrid[i];
										Qabs += dQabs;
										if(dQabs>1e-12*Qabs) { isNonzero=true; iLastNZ=i+2; }
									}
									if(!isNonzero) continue;
									if(iLastNZ < nGrid)
									{	Qijl.r.resize(iLastNZ);
										Qijl.dr.resize(iLastNZ);
										Qijl.f.resize(iLastNZ);
									}
									//Store in Qradial:
									QijIndex qIndex = { l1, p1, l2, p2, l };
									Qijl.transform(l, dG, nGridLoc, Qradial[qIndex]);
									//Store Qint = integral(Qradial) when relevant:
									if(l1==l2 && !l)
									{	double Qint_ij = Qijl.transform(0,0)/(4*M_PI);
										Qint[l1].data()[Qint[l1].index(p1,p2)] = Qint_ij;
										Qint[l1].data()[Qint[l1].index(p2,p1)] = Qint_ij;
									}
								}
							}
						}
					}
					if(!Qradial.size()) Qint.clear(); //Norm-conserving pseudopotentials that had been written in Ultrasoft format for some reason
				}
			}
		}
		else if(tag.name == "PP_PSWFC")
		{	lPsi.assign(nPsi, -1); //angular momentum per atomic orbital
			std::vector<RadialFunctionR> psi(nPsi, RadialFunctionR(nGrid)); //orbitals
			std::vector<double> eigs(nPsi); //eigenvalues
			bool haveEigs = true;
			while(true)
			{	XMLtag tagPsi(is);
				if(!tagPsi.name.length()) break;
				if(tagPsi.name.substr(0,7) == "PP_CHI.") //Atomic orbital
				{	int iPsi = atoi(tagPsi.name.substr(7).c_str());
					if(iPsi<1 || iPsi>nPsi)
						die("  Invalid orbital index %d (not in range [1,%d]) in section '%s'\n", iPsi, nPsi, tagPsi.name.c_str());
					iPsi--; //to zero-based index
					//Angular momentum
					lPsi[iPsi] = atoi(tagPsi.getAttribute("l").c_str());
					if(lPsi[iPsi]<0)
						die("  Invalid projector angular momentum %d (not >= 0) in section '%s'\n", lPsi[iPsi], tagPsi.name.c_str());
					//Eigenvalue and general info:
					double occ = atof(tagPsi.getAttribute("occupation").c_str());
					string label = tagPsi.getAttribute("label");
					eigs[iPsi] = atof(tagPsi.getAttribute("pseudo_energy", false).c_str()) * 0.5; //convert from Ry to Eh
					logPrintf("    %-3s   l: %d   occupation: %4.1lf", label.c_str(), lPsi[iPsi], occ);
					if(eigs[iPsi])
						logPrintf("   eigenvalue: %lf", eigs[iPsi]);
					else
						haveEigs=false;
					logPrintf("\n");
					//Orbital function:
					psi[iPsi].set(rGrid, drGrid);
					psi[iPsi].f = tagPsi.readData(nGrid);
				}
				else
				{	logPrintf("  NOTE: ignored section '%s'\n", tagPsi.name.c_str());
					tagPsi.ignoreAndClose();
				}
				tagPsi.close();
			}
			if(nPsi > 0)
			{	logPrintf("  Transforming atomic orbitals to a uniform radial grid of dG=%lg with %d points.\n", dG, nGridNL);
				psiRadial.resize(lMax+1);
				if(haveEigs) atomEigs.resize(lMax+1);
				for(int iPsi=0; iPsi<nPsi; iPsi++)
				{	int l = lPsi[iPsi];
					if(l>lMax)
					{	atomEigs.resize(l+1);
						psiRadial.resize(l+1);
					}
					for(int i=0; i<nGrid; i++)
						psi[iPsi].f[i] *= (rGrid[i] ? 1./rGrid[i] : 0);
					psiRadial[l].push_back(RadialFunctionG());
					psi[iPsi].transform(l, dG, nGridNL, psiRadial[l].back());
					if(haveEigs) atomEigs[l].push_back(eigs[iPsi]);
				}
			}
		}
		else if(tag.name == "PP_SPIN_ORB")
		{	jNL.resize(nBeta, -1.);
			jPsi.resize(nPsi, -1.);
			while(true)
			{	XMLtag tagSO(is);
				if(!tagSO.name.length()) break;
				if(tagSO.name.substr(0,11) == "PP_RELBETA.") //Projector SO info
				{	int iBeta = atoi(tagSO.name.substr(11).c_str());
					if(iBeta<1 || iBeta>nBeta)
						die("  Invalid projector index %d (not in range [1,%d]) in section '%s'\n", iBeta, nBeta, tagSO.name.c_str());
					iBeta--; //to zero-based index
					jNL[iBeta] = atof(tagSO.getAttribute("jjj").c_str());
				}
				else if(tagSO.name.substr(0,10) == "PP_RELWFC.") //Atomic orbital SO info
				{	int iPsi = atoi(tagSO.name.substr(10).c_str());
					if(iPsi<1 || iPsi>nPsi)
						die("  Invalid orbital index %d (not in range [1,%d]) in section '%s'\n", iPsi, nPsi, tagSO.name.c_str());
					iPsi--; //to zero-based index
					jPsi[iPsi] = atof(tagSO.getAttribute("jchi").c_str());
				}
				else
				{	logPrintf("  NOTE: ignored section '%s'\n", tagSO.name.c_str());
					tagSO.ignoreAndClose();
				}
				tagSO.close();
			}
		}
		else //Unused/unknown section (ignore, but mention in log file just in case)
		{	logPrintf("  NOTE: ignored section '%s'\n", tag.name.c_str());
			tag.ignoreAndClose();
		}
		tag.close();
	}
	tagUPF.ignoreAndClose();
	
	//Process j's for relativistic pseudopotentials
	if(jNL.size())
	{	Vnl2j.resize(VnlRadial.size());
		psi2j.resize(psiRadial.size());
		for(int iBeta=0; iBeta<nBeta; iBeta++)
		{	int j2 = round(2*jNL[iBeta]), l = lNL[iBeta];
			if(not (j2==2*l-1 || j2==2*l+1))
				die("  Total angular momentum %lg incompatible with orbital angular momentum %d for projector# %d\n", jNL[iBeta], l, iBeta+1);
			Vnl2j[l].push_back(j2);
		}
		for(int iPsi=0; iPsi<nPsi; iPsi++)
		{	int j2 = round(2*jPsi[iPsi]), l = lPsi[iPsi];
			if(not (j2==2*l-1 || j2==2*l+1))
				die("  Total angular momentum %lg incompatible with orbital angular momentum %d for orbital# %d\n", jPsi[iPsi], l, iPsi+1);
			psi2j[l].push_back(j2);
		}
	}
	
	coreRadius = *coreRadii.rbegin(); //max of all core radii above (used for overlap checks during geometry opt)
}
