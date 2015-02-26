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

#include <commands/command.h>
#include <electronic/Everything.h>
#include <core/string.h>
#include <limits.h>
#include <algorithm>

using std::map;
using std::pair;
using std::make_pair;
using std::set;

//! Process the statically initialized command map, add dependency info
//! and check for inconsistencies (static ones i.e. bugs in the code, not input file)
struct ProcessedCommandMap : public map<string, pair<int,Command*> >
{
	int nPasses; //!< maximum number of passes required based on dependency tree

	//! Calls getCommandMap(), checks and encodes dependency information
	//! commands paired with integer 0 have no dependencies,
	//! those with 1 depend only on others with 0 etc.
	ProcessedCommandMap()
	{	const map<string,Command*>& origMap = getCommandMap();
		//Copy into *this with all the integers set to INT_MAX
		for(map<string,Command*>::const_iterator i=origMap.begin(); i!=origMap.end(); i++)
			(*this)[i->first] = make_pair(INT_MAX, i->second);
		//Make multiple passes till all the commands have been handled:
		for(nPasses=0; ; nPasses++)
		{	int nProcessed = 0; //number processed in this pass
			int nRemaining = 0; //number not yet processed
			for(iterator i=begin(); i!=end(); i++)
				if(i->second.first>nPasses) //not yet processed
				{	Command& ci = *(i->second.second);
					bool ready = true;
					//Loop through all its required commands:
					for(set<string>::iterator r=ci.requires.begin(); r!=ci.requires.end(); r++)
					{	iterator j = find(*r);
						if(j==end())
							die("Inconsistency in command map: '%s' in '%s's requires list is not a valid command.\n",
								r->c_str(), ci.name.c_str());
						if(j->second.first>=nPasses) { ready=false; break; } //unsatisfied dependency
					}
					if(ready)
					{	i->second.first = nPasses;
						nProcessed++;
					}
					else nRemaining++;
				}
			if(!nRemaining) break; //all commands have been processed
			if(!nProcessed)
			{	logPrintf("There is a cyclic dependency somewhere amongst the following commands:\n");
				for(iterator i=begin(); i!=end(); i++)
					if(i->second.first>nPasses)
						logPrintf("\t%s\n", i->second.second->name.c_str());
				die("Cyclic dependency in command map.\n");
			}
		}
		//Make forbids mutual:
		typedef map<string,Command*>::const_iterator Citer;
		for(Citer i=origMap.begin(); i!=origMap.end(); i++)
		{	Command& ci = *(i->second);
			for(set<string>::iterator f=ci.forbids.begin(); f!=ci.forbids.end(); f++)
			{	Citer j = origMap.find(*f);
				if(j==origMap.end())
					die("Inconsistency in command map: '%s' in '%s's forbids list is not a valid command.\n",
						f->c_str(), ci.name.c_str());
				Command& cj = *(j->second);
				if(!cj.forbids.count(ci.name) && cj.section==ci.section) //don't report cross section errors as they can't be fixed
				{	logPrintf("Command map WARNING: '%s' forbids '%s', but not vice versa (making this so)\n",
						ci.name.c_str(), cj.name.c_str());
					cj.forbids.insert(ci.name);
				}
			}
		}
		//Check for inconsistency in the forbids:
		for(Citer i=origMap.begin(); i!=origMap.end(); i++)
		{	Command& ci = *(i->second);
			set<Command*> deps, forbids;
			getDependencies(ci, origMap, deps, forbids);
			//Find the intersection of the dependencies and the forbidden commands:
			std::vector<Command*> common(std::max(deps.size(),forbids.size()));
			common.erase(
				set_intersection(deps.begin(),deps.end(), forbids.begin(),forbids.end(), common.begin()),
				common.end() );
			if(common.size())
			{	logPrintf("Command '%s' (in)directly requires as well as forbids the following commands:\n", ci.name.c_str());
				for(unsigned j=0; j<common.size(); j++)
					logPrintf("\t%s\n", common[j]->name.c_str());
				die("Forbid/require inconsistency in command map.\n");
			}
		}
	}

private:
	//Add all the dependencies of command c in dep,
	//and all the commands forbidden by c and its dependencies to forbids
	void getDependencies(const Command& c, const map<string,Command*>& cmdMap, set<Command*>& dep, set<Command*>& forbids)
	{	for(set<string>::iterator r=c.requires.begin(); r!=c.requires.end(); r++)
		{	Command* cr = cmdMap.find(*r)->second;
			dep.insert(cr);
			for(set<string>::iterator f=cr->forbids.begin(); f!=cr->forbids.end(); f++)
				forbids.insert(cmdMap.find(*f)->second);
			getDependencies(*cr, cmdMap, dep, forbids);
		}
	}
};

//! Print buf to logPrintf, but make sure every line is commented:
void printCommented(string buf)
{	istringstream iss(buf);
	while(!iss.eof())
	{	//read a line
		string line; getline(iss, line);
		//Check if the first non white-space char is #:
		char firstChar = ' ';
		sscanf(line.c_str(), "%c", &firstChar);
		if(firstChar!='#') logPrintf("# "); //add the comment character
		logPrintf("%s\n", line.c_str());
	}
}

//! Replace all ${...} sequences in s with values of corresponding environment variables
void environmentSubstitute(string& s)
{	while(true)
	{	size_t varStart = s.find("${");
		if(varStart==string::npos) break; //no variables to replace
		size_t varStop = s.find("}", varStart);
		if(varStop==string::npos)
			die("Unterminated environment variable '%s'.\n", s.substr(varStart, varStop-varStart).c_str());
		string varName = s.substr(varStart+2,varStop-(varStart+2));
		const char* varValue = getenv(varName.c_str()); //returns null pointer if not found
		s.replace(varStart, varStop+1-varStart, varValue ? varValue : "");
	}
}

//! Read line from stream, after ignoring comments and handling line continuation
string readLine(istream& is)
{	string line;
	bool continueLine = true;
	while(continueLine && !is.eof())
	{	string lineIn; getline(is, lineIn);
		//Get rid of comments
		size_t commentStart = lineIn.find('#');
		if(commentStart != string::npos) //found a #
			lineIn.erase(commentStart); //remove everything from the # onwards
		//Remove leading and trailing whitespace:
		trim(lineIn);
		environmentSubstitute(lineIn);
		//Check for line continuation:
		continueLine = false;
		if(lineIn.length() && lineIn[lineIn.length()-1]=='\\')
		{	continueLine = true;
			lineIn.erase(lineIn.length()-1);
		}
		//Add to the line buffer:
		if(line.length()) line += ' '; //Add a space between joined lines
		line += lineIn;
	}
	return line;
}

//! Read file filename.back() and add pairs of commands and parameters to input
//! Handle line continuation characters and discard comments in the process.
//! Also process include definitions, and use the contents of filename to detect cyclic include
void readInputFile(std::vector<string>& filename, std::vector< pair<string,string> >& input)
{	//Open last file in list (or stdin if its blank)
	ifstream ifs;
	if(filename.back().length())
	{	ifs.open(filename.back().c_str()); //open the last file in the list
		if(!ifs.is_open()) die("Could not open file '%s' for reading.\n", filename.back().c_str());
	}
	else
		logPrintf("Waiting for commands from stdin (end input with EOF (Ctrl+D)):\n");
	istream& is = filename.back().length() ? ifs : std::cin;
	//Read the file line by line:
	while(!is.eof())
	{	string line = readLine(is);
		if(!line.length()) continue; //ignore empty lines
		//Break the line into command name and parameter list:
		string cmd = line.substr(0, line.find_first_of(" \t\n\r"));
		string params = line.substr(cmd.length());
		//Handle includes:
		if(cmd=="include")
		{	trim(params);
			if(!params.length()) die("No file-name specified for command include.\n");
			//Check for cyclic includes:
			if(find(filename.begin(), filename.end(), params) != filename.end())
				die("File '%s' tries to (in)directly include itself.\n", params.c_str());
			//Read included file recursively:
			filename.push_back(params);
			readInputFile(filename, input);
		}
		else input.push_back(make_pair(cmd,params));
	}
	//Done with this file, remove it from the list (stack, really) of files
	filename.pop_back();
}

std::vector< pair<string,string> > readInputFile(string filename)
{	std::vector< pair<string,string> > input;
	//Read input on head:
	if(mpiUtil->isHead())
	{	std::vector<string> filenameList(1, filename);
		readInputFile(filenameList, input);
	}
	//Broadcast to other processes
	if(mpiUtil->nProcesses()>1)
	{	//Synchronize lengths:
		int nInputs = input.size();
		mpiUtil->bcast(nInputs);
		if(!mpiUtil->isHead()) input.resize(nInputs);
		//Serialize to string and synchronize content:
		string inputStr;
		if(mpiUtil->isHead())
		{	ostringstream oss;
			for(const auto& p: input)
				oss << p.first << '\n' << p.second << '\n';
			inputStr = oss.str();
		}
		mpiUtil->bcast(inputStr);
		if(!mpiUtil->isHead())
		{	istringstream iss(inputStr);
			for(auto& p: input)
			{	getline(iss, p.first);
				getline(iss, p.second);
			}
		}
	}
	return input;
}

//! Call Command::process with error handling, count updating etc:
void safeProcess(Command& c, string params, Everything& everything,
	map<string,int>& encountered, std::vector< pair<Command*,string> >& errors)
{
	//Update the encountered count:
	encountered[c.name]++;
	//Check dependencies:
	string unsatisfiedDeps;
	for(set<string>::iterator r=c.requires.begin(); r!=c.requires.end(); r++)
		if(!encountered[*r])
			unsatisfiedDeps += ' '+ (*r);
	if(unsatisfiedDeps.length())
	{	errors.push_back(make_pair(&c, "is missing dependencies {"+unsatisfiedDeps+" }") );
		return;
	}
	//Check for multiplicity:
	if(!c.allowMultiple && encountered[c.name]>1)
	{	errors.push_back(make_pair(&c, "should not be issued more than once") );
		return;
	}
	//Check forbidden list:
	for(const string& forbidden: c.forbids)
		if(encountered[forbidden])
			errors.push_back(make_pair(&c, " is incompatible with command " + forbidden) );
	
	//Try running the command
	try
	{	ParamList pl(params);
		c.process(pl, everything);
		//Check remainder:
		string remainder = pl.getRemainder();
		if(remainder.length())
			throw string("Extra arguments '" + remainder + "'  at end of command");
	}
	catch(string err)
	{	errors.push_back(make_pair(&c, " with command line:\n"
			"\t"+c.name+params+"\n"
			"failed with message:\n\t"+err+"\n"
			"The correct usage of this command is described below:\n"
			"\t"+c.name+' '+c.format+'\n'+c.comments ) );
		return;
	}
}

void parse(std::vector< pair<string,string> > input, Everything& everything, bool printDefaults)
{
	//Process the input in multiple passes: (Be helpful by collecting all errors, if any, before dying)
	ProcessedCommandMap cmap;
	set<string> unknown; //unknown command names
	map<string,int> encountered; //command names, and the number of times they were encountered
	std::vector< pair<Command*,string> > errors; //command classes that encountered errors, with error messages
	
	//First check for, take note of, and remove unknown commands
	for(unsigned i=0; i<input.size();)
	{	if(cmap.find(input[i].first)==cmap.end())
		{	unknown.insert(input[i].first);
			input.erase(input.begin()+i);
		}
		else i++;
	}
	
	//All command names in input are now known to be in cmap (so can safely use cmap[] instead of cmap.find())
	for(int pass=0; pass<=cmap.nPasses; pass++)
	{	//Run through all commands marked for this pass:
		for(unsigned i=0; i<input.size(); i++)
		{	pair<int,Command*>& icPair = cmap[input[i].first];
			if(icPair.first==pass)
			{	Command& c = *(icPair.second);
				//Check for empty parameter list:
				string trimmedParams = input[i].second;
				trim(trimmedParams);
				if(!trimmedParams.length() && c.emptyParamError.length())
					errors.push_back(make_pair(&c,
						" cannot be issued without any parameters:\n"
						+ c.emptyParamError
						+ (c.hasDefault ? "\nOmit this command altogether for the default behaviour.\n" : "")));
				//Process command:
				string params = input[i].second + " "; //add space at end to prevent eof on last parameter
				safeProcess(c, params, everything, encountered, errors);
			}
		}
		//Run defaults where necessary:
		for(ProcessedCommandMap::iterator i=cmap.begin(); i!=cmap.end(); i++)
			if(i->second.first==pass)
			{	Command& c = *(i->second.second);
				if(c.hasDefault && !encountered[c.name])
					safeProcess(c, "", everything, encountered, errors);
			}
	}
	//Quit on errors:
	size_t errTot = unknown.size() + errors.size();
	if(errTot)
	{	logPrintf("\n\nEncountered %lu errors while parsing input:\n\n", errTot);
		for(set<string>::iterator i=unknown.begin(); i!=unknown.end(); i++)
			logPrintf("'%s' is not a valid command.\n", i->c_str());
		logPrintf("\n");
		for(unsigned i=0; i<errors.size(); i++)
			logPrintf("Command '%s' %s\n\n", errors[i].first->name.c_str(), errors[i].second.c_str());
		die("\n\nInput parsing failed with %lu errors (run with -t for command syntax)\n\n", errTot);
	}
	//Print status:
	std::map<string,int> explicitlyEncountered; //List of commands explicitly in input file (with corresponding multiplicities)
	if(!printDefaults) for(auto cmd: input) explicitlyEncountered[cmd.first]++;
	logPrintf("\n\nInput parsed successfully to the following command list (%sincluding defaults):\n\n", printDefaults ? "" : "not ");
	for(auto i: (printDefaults ? encountered : explicitlyEncountered))
	{	for(int iRep=0; iRep<i.second; iRep++) //handle repetitions
		{	Command& c = *(cmap[i.first].second);
			logPrintf("%s ", c.name.c_str());
			c.printStatus(everything, iRep);
			logPrintf("\n");
		}
	}
	logPrintf("\n\n");
}

inline void processDefaults(Everything& everything, ProcessedCommandMap& cmap)
{	//Run through the commands that have defaults in dependency order
	for(int pass=0; pass<=cmap.nPasses; pass++)
		for(ProcessedCommandMap::iterator i=cmap.begin(); i!=cmap.end(); i++)
			if(i->second.first==pass)
			{	Command& ci = *(i->second.second);
				if(ci.hasDefault)
				{	try
					{	ParamList pl(""); //run with no arguments
						ci.process(pl, everything);
					}
					catch(string str)
					{	die("BUG in command '%s'; default call raised error:\n\t%s\n",
							ci.name.c_str(), str.c_str());
					}
				}
			}
}

inline string describeSyntax()
{	return
		" * The input file may contain commands in any order; commands will be\n"
		"   automatically processed in an order that satisfies all dependencies\n"
		"\n"
		" * Each command is a single line, but lines can be continued using '\\'\n"
		"\n"
		" * Whitespace is unimportant, except for separating words in each line\n"
		"\n"
		" * Any text following '#' on an input line is treated as a comment\n"
		"\n"
		" * 'include <file>' can be used to include commands from <file>\n"
		"\n"
		" * Each instance of ${xxx} is replaced by environment variable 'xxx'\n"
		"   (Variable substitution occurs before command/include processing)\n";
}

void printDefaultTemplate(Everything& everything)
{	ProcessedCommandMap cmap;
	processDefaults(everything, cmap);
	//Now print the comments and status of all commands in map (alphabetical) order:
	for(ProcessedCommandMap::iterator i=cmap.begin(); i!=cmap.end(); i++)
	{	Command& ci = *(i->second.second);
		//Print format and comments:
		printCommented(ci.name+' '+ci.format+'\n'+ci.comments);
		//Print status:
		if(ci.hasDefault)
		{	logPrintf("%s ", ci.name.c_str());
			ci.printStatus(everything, 0);
			logPrintf("\n");
		}
		logPrintf("\n");
	}
	//Print general input file help:
	logPrintf("\n\n"
		"# +------------------------------------------------------------------------+\n"
		"# |                       JDFTx input file format                          |\n"
		"# +------------------------------------------------------------------------+\n"
		"#\n"
	);
	printCommented(describeSyntax());
}

inline string htmlEscapeCharacters(string s)
{	string sOut; sOut.reserve(s.length()+100);
	for(const char& c: s)
	{	switch(c)
		{	case '<': sOut.append("&lt;"); break;
			case '>': sOut.append("&gt;"); break;
			default: sOut.push_back(c); break;
		}
	}
	return sOut;
}

inline string commandNameToID(string name)
{	string id = "Command";
	bool prevCaps = false, curCaps = true;
	for(const char c: name)
	{	if(c == '-')
		{	curCaps = !prevCaps;
		}
		else
		{	char cOut = curCaps ? toupper(c) : c;
			id += cOut;
			curCaps = false;
			prevCaps = isupper(cOut);
		}
	}
	return id;
}


inline string htmlAddLinks(string s)
{	static std::set<string> cExcluded;
	if(!cExcluded.size()) //Add commands to exclude from auto-linking (most likely because they are common words)
	{	cExcluded.insert("lattice"); cExcluded.insert("ion"); cExcluded.insert("basis");
		cExcluded.insert("fluid"); cExcluded.insert("wavefunction"); cExcluded.insert("symmetries");
		cExcluded.insert("debug"); cExcluded.insert("dump"); cExcluded.insert("polarizability");
		cExcluded.insert("vibrations"); 
	}
	std::map<string,Command*> cmap = getCommandMap();
	string sOut; sOut.reserve(s.length()+400);
	const char* delim = " \t\n.,;:)([]!?'\""; //delimiters for tokenization
	size_t pos = 0;
	string wordPrev, wordPrev2;
	while(pos<s.length())
	{	//Forward all delimiters unmodified:
		size_t posNext = std::min(s.find_first_not_of(delim, pos), s.length());
		sOut.append(s.substr(pos,posNext-pos));
		pos = posNext;
		//Check next word and add link if it is a command name:
		posNext = std::min(s.find_first_of(delim, pos), s.length());
		string word = s.substr(pos,posNext-pos);
		if(!(wordPrev2=="see" && wordPrev=="command")
			&& (cmap.find(word)==cmap.end() || cExcluded.count(word))) //Not a command name or excluded command, and not explicitly linked
			sOut.append(word); //just forward
		else //Command name that is not to be excluded
			sOut.append("\\ref " + commandNameToID(word) + " \"" + word + "\""); //add a link
		pos = posNext;
		wordPrev2 = wordPrev;
		wordPrev = word;
	}
	return sOut;
}

inline void printHTMLformatted(string s)
{	fputs(htmlAddLinks(htmlEscapeCharacters(s)).c_str(), globalLog);
}

void writeCommandManual(Everything& everything, string section)
{	ProcessedCommandMap cmap;
	processDefaults(everything, cmap);
	//Print a doxygen page for all commands:
	logPrintf("//Auto-generated using writeCommandManual(everything, \"%s\")\n", section.c_str());
	logPrintf("//Do not edit manually: instead edit the documentation strings in the code.\n");
	for(auto& i: cmap)
	{	Command& ci = *(i.second.second);
		if(ci.section != section) continue; //not in current section/executable
		//Start page:
		string id = commandNameToID(ci.name);
		logPrintf("\n/** \\page %s %s\n", id.c_str(), ci.name.c_str());
		//Print syntax:
		logPrintf("Syntax:\n");
		logPrintf("-------\n");
		logPrintf("\n    %s %s\n", ci.name.c_str(), ci.format.c_str());
		logPrintf("\n");
		if(ci.section != "jdftx")
			logPrintf("<b>Note:</b> only available in calculations using the '%s' executable.\n", ci.section.c_str());
		//Print description (comments):
		logPrintf("Description:\n");
		logPrintf("------------\n");
		printHTMLformatted(ci.comments+'\n');
		logPrintf("\n");
		//Print properties:
		logPrintf("Properties:\n");
		logPrintf("-----------\n");
		//--- requires
		logPrintf("\n<b>Requires:</b>\n");
		for(const string& s: ci.requires)
			logPrintf("\\bigsep \\ref %s \"%s\"\n", commandNameToID(s).c_str(), s.c_str());
		if(!ci.requires.size()) logPrintf("\\bigsep (None)\n");
		//--- forbids
		logPrintf("\n<b>Forbids:</b>\n");
		for(const string& s: ci.forbids)
			logPrintf("\\bigsep \\ref %s \"%s\"\n", commandNameToID(s).c_str(), s.c_str());
		if(!ci.forbids.size()) logPrintf("\\bigsep (None)\n");
		//--- allowMultiple
		logPrintf("\n<b>Allow multiple:</b>\\bigsep %s\n", boolMap.getString(ci.allowMultiple));
		//--- default
		logPrintf("\n<b>Default:</b>\n");
		if(ci.hasDefault)
		{	logPrintf("\n    %s ", ci.name.c_str());
			ci.printStatus(everything, 0);
			logPrintf("\n");
		}
		else logPrintf("\\bigsep (None)\n");
		//Link back to main page and close:
		logPrintf("\nBack to: \\ref Commands or \\ref CommandIndex \n");
		logPrintf("*/\n");
	}
}
