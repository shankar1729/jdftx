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

#ifndef JDFTX_COMMAND_COMMAND_H
#define JDFTX_COMMAND_COMMAND_H

#include <commands/ParamList.h>
#include <core/string.h>
#include <memory>
#include <set>

//! @addtogroup Input
//! @{

/** @file command.h
@brief Provides the base class and various helpers for defining commands in the input file 

To create a new command, paste this code block into a cpp file,
replace all occurences of dummy with something suitable, add the
documentation strings (format and comment) and the parse code.
See any of the existing commands for examples of using ParamList effectively.
Note that there is no need to edit any global list of commands, this happens
automatically during static initialization, just create an object of the
new derived Command class at file scope. However, for coherent documentation,
edit doc/commands.dox and link to the command in the appropriate section.

	struct CommandDummy : public Command
	{
		CommandDummy() : Command("dummy")
		{
			format = "";
			comments = "";
		}

		void process(ParamList& pl, Everything& e)
		{
		}

		void printStatus(Everything& e, int iRep)
		{
		}
	}
	commandDummy;
*/

//! @brief Abstract base class for all commands
class Command
{
public:
	string name; //!< Identifier for the command in the input file. Must be unique!
	string format; //!< Usage syntax for the command (excluding the command name)
	string comments; //!< Detailed help for the command which goes into jdftx -t as well as the Doxygen manual. Please check formatting in both versions.
	string section; //!< Which executable the command belongs to, and hence which section it must be documented under.
	string category; //!< Category of command under which to list the documentation
	string subcategory; //!< Subcategory of command under which to list the documentation
	
	std::set<string> requires; //!< Names of other commands that this one requires; those commands will be processed before this one.
	std::set<string> forbids; //!< Names of other commands that this one is incompatibile with.

	bool allowMultiple; //!< Whether this command can occur multiple times in an input file (default: false).
	bool hasDefault; //!< Whether this command has a default (default: false). If true, process() will be called with an empty parameter list.
	string emptyParamError; //!< Error message if command is manually issued with no parameters. Default: empty i.e. allow empty parameter list.
	
	//! @brief Process the command from its command line
	//! @param pl Parser utility object that simplifies extraction of parameters from the commandline
	//! @param e Object of class Everything into which the settings extracted form this command should be stored
	//! @throw string On any syntax / parse errors, throw a string containing a brief error (preferably single phrase) error message.
	virtual void process(ParamList& pl, Everything& e)=0;

	//! @brief Print a command line that would result in the current status
	//! @param e Object of class Everything whose status to reproduce
	//! @param iRep For commands with allowMultiple=true, printStatus will be called multiple times and iRep is the index of the present call
	virtual void printStatus(Everything& e, int iRep)=0;

protected:
	//! This base class constructor adds the current command to a map from names to Command pointers
	//! which can be accessed using getCommandMap(). This enables safe static initialization of the command list.
	//! @param name Unique name of the command
	//! @param path Documentation path for the command in the format
	//!   "section/category/subcategory", where section is the name
	//!   of the executable (jdftx, phonon or wannier),
	//!   and the remaining help organize commands in the doc
	Command(string name, string path);

	void require(string); //!< utility to add a command to the requires list
	void forbid(string); //!< utility to add a command to the forbids list
};

//! Base class for a deprecated command which will translate old syntax into the new command that replaces it
class DeprecatedCommand
{
public:
	string name;
	DeprecatedCommand(string name);
	
	//Return a replacement command and arguments pair
	virtual std::pair<string,string> replace(ParamList& pl) const=0;
};

std::map<string,Command*>& getCommandMap(); //!< Retrieve the map from command names to objects created during static initialization
std::map<string,DeprecatedCommand*>& getDeprecatedMap(); //!< Retrieve the map from deprecated command names to objects

static EnumStringMap<bool> boolMap(false, "no", true, "yes"); //!< utility to parse yes/no into bools


//! @brief Process the EnumStringMap::optionList() to add descriptions using an arbitrary functor
//! @tparam GetDescription Function/functor with the signature `string GetDescription(const string&)`.
//! @param optionList List of strings separated by pipe characters as returned by EnumStringMap::optionList()
//! @param getDescription A function/functor that returns the description given an option from optionList
//! @param spacer Spacer to insert before each option/description pair
//! @return string containing the options along with their descriptions
template<typename GetDescription>
string addDescriptions(string optionList, const GetDescription& getDescription, string spacer="\n+ ")
{
	//Determine max width of name, so as to align:
	istringstream iss(optionList);
	size_t nameWidth=0;
	while(!iss.eof())
	{	string name;
		getline(iss, name, '|');
		trim(name);
		nameWidth = std::max(name.length(),nameWidth);
	}
	//Process decsription list:
	iss.seekg(0, std::ios::beg); //rewind;
	string ret;
	while(!iss.eof())
	{	//Get a name from the option list:
		string name;
		getline(iss, name, '|');
		trim(name);
		if(!name.length()) break;
		//Get the description:
		string desc = getDescription(name);
		//Pad the name to required width:
		if(name.length()<nameWidth)
			name.resize(nameWidth, ' ');
		ret += spacer + name;
		if(desc.length()) ret += ": " + desc;
	}
	return ret;
}

//! @brief Null description function that can be used with addDescriptions().
//! This can be useful for simply reformatting EnumStringMap::optionList() 
inline string nullDescription(const string&)
{	return string();
}

//! @brief Generate a description functor for addDescriptions() from an EnumStringMap
//! @tparam Enum The template parameter of EnumStringMap
template<typename Enum> struct LinkDescription
{	const EnumStringMap<Enum>& nameMap;
	const EnumStringMap<Enum>& descMap;
	
	//! @param nameMap Link between enum members to corresponding name strings
	//! @param descMap Link between enum members to corresponding description strings
	LinkDescription(const EnumStringMap<Enum>& nameMap, const EnumStringMap<Enum>& descMap)
	: nameMap(nameMap), descMap(descMap)
	{
	}
	
	//! @brief Look name up in nameMap, and return the corresponding description in descMap
	//! @param name One of the string entries in nameMap
	//! @return The corresponding string in descMap
	string operator()(const string& name) const
	{	Enum type = Enum();
		bool nameFound = nameMap.getEnum(name.c_str(), type);
		assert(nameFound);
		return descMap.getString(type);
	}
};

//! @brief Helper function to select appropriate LinkDescription<> by overloading
//! @tparam Enum The template parameter of EnumStringMap
//! @param nameMap Link between enum members to corresponding name strings
//! @param descMap Link between enum members to corresponding description strings
template<typename Enum>
LinkDescription<Enum> linkDescription(const EnumStringMap<Enum>& nameMap, const EnumStringMap<Enum>& descMap)
{	return LinkDescription<Enum>(nameMap, descMap);
}

//! @brief Find species matching an id (and create it from a wildcard if necessary)
//! @param id Psuedopotential species identifier string (typically chemical symbol)
//! @param e Reference to Everything
//! @return Smart pointer to SpeciesInfo corresponding to id.
std::shared_ptr<class SpeciesInfo> findSpecies(string id, Everything& e);

//! @brief Check if file is readable in an MPI friendly way.
//! This function must be called from all processes,
//! the file will be checked on one process 
//! and the result broadcast to all processes.
//! @param fname File name to check for
//! @return Whether file is readable
bool isReadable(string fname);

//! Set target to the filename for reading variable varName from file with pattern
//! specified by filenamePattern if that file exists and is readable
void setAvailableFilename(string filenamePattern, string varName, string& target);

//! Apply setAvailableFilename for all standard input variables (action of CommandInitialState)
void setAvailableFilenames(string filenamePattern, Everything& e);

//! @}
#endif // JDFTX_COMMAND_COMMAND_H
