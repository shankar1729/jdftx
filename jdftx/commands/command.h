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

#include <electronic/common.h>
#include <commands/ParamList.h>
#include <core/string.h>
#include <set>

//!@file command.h
/**  To create a new command, paste this code block into a cpp file,
replace all occurences of dummy with something suitable, add the
documentation strings (format and comment) and the parse code.
See any of the existing commands for examples of using ParamList effectively.
Note that there is no need to edit any global list of commands, this happens
automatically during static initialization, just create an object of the
new derived Command class at file scope.

//! Implements command "dummy"
typedef struct  CommandDummy : public Command
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

//! Abstract base class for all commands
class Command
{
public:
	string name; //!< identifier for the command in the input file (MUST be unique!)
	string format; //!< usage syntax for the command (excluding the command name)
	string comments; //!< more detailed help for the command (note comment characters are automatically added!)

	std::set<string> requires; //!< names of other commands this requires (those will be processed before)
	std::set<string> forbids; //!< names of other commands this is incompatibile with

	bool allowMultiple; //!< whether this command can occur multiple times in an input file (default: false)
	bool hasDefault; //!< whether this command has a default (process is called with no parameters) (default: false)
	bool overallAction; //!< whether this command requires an overall action to be performed (default: false)
	string emptyParamError; //!< error message if command is manually issued with no parameters (default: "" => allow empty parameter list)
	
	//! Process the command given its parameters.
	//! Throw a string with the message as an exception on error
	virtual void process(ParamList& pl, Everything& e)=0;

	//! Print the parameter set that would have reuslted in the current status
	//! For commands with allowMultiple=true, this will be called with varying iRep
	//! If the command takes no parameters, then nothing should be printed
	virtual void printStatus(Everything& e, int iRep)=0;

protected:
	//! This base class constructor adds the current command to a map from names to Command pointers
	//! which can be accessed using getCommandMap(). This enables safe static initialization of the command list.
	Command(string name);

	void require(string); //!< utility to add a command to the requires list
	void forbid(string); //!< utility to add a command to the forbids list
};

std::map<string,Command*>& getCommandMap(); //!< Get the list of commands created by static initialization

static EnumStringMap<bool> boolMap(false, "no", true, "yes"); //!< utility to parse yes/no into bools

//! Process the EnumStringMap::optionList() to add descriptions using an arbitrary functor
template<typename GetDescription>
string addDescriptions(string optionList, const GetDescription& getDescription, string spacer="\n   ")
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
		ret += spacer + name + " " + desc;
	}
	return ret;
}

//! Generate a description functor using a EnumStringMap holding the decsriptions:
template<typename Enum> struct LinkDescription
{	const EnumStringMap<Enum>& nameMap;
	const EnumStringMap<Enum>& descMap;
	
	LinkDescription(const EnumStringMap<Enum>& nameMap, const EnumStringMap<Enum>& descMap)
	: nameMap(nameMap), descMap(descMap)
	{
	}
	
	//! Look name up in nameMap, and return the corresponding description in descMap
	string operator()(const string& name) const
	{	Enum type = Enum();
		assert(nameMap.getEnum(name.c_str(), type));
		return descMap.getString(type);
	}
};

//! Helper function to select appropriate LinkDescription<> by overloading
template<typename Enum>
LinkDescription<Enum> linkDescription(const EnumStringMap<Enum>& nameMap, const EnumStringMap<Enum>& descMap)
{	return LinkDescription<Enum>(nameMap, descMap);
}

#endif // JDFTX_COMMAND_COMMAND_H
