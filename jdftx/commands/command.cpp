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

using std::map;

map<string,Command*>& updateCommandMap(Command* command=0)
{	static map<string,Command*> commandMap; //make this local static for correct initialization order
	if(command) commandMap[command->name] = command;
	return commandMap;
}
map<string,Command*>& getCommandMap()
{	return updateCommandMap(0);
}

void fixCategoryName(string& name)
{	//Replace spaces with underscore:
	for(char& c: name) if(c==' ') c='_';
	//Replace empty string with NULL:
	if(!name.length()) name = "NULL";
}

Command::Command(string name, string path)
: name(name), allowMultiple(false), hasDefault(false)
{	updateCommandMap(this);
	
	//Parse the path to get section, category and subcategory:
	istringstream iss(path);
	getline(iss, section, '/');
	getline(iss, category, '/'); fixCategoryName(category);
	getline(iss, subcategory, '/'); fixCategoryName(subcategory);
}

void Command::require(string name) { requires.insert(name); }
void Command::forbid(string name) { forbids.insert(name); }


map<string,DeprecatedCommand*>& updateDeprecatedMap(DeprecatedCommand* dc=0)
{	static map<string,DeprecatedCommand*> deprecatedMap; //make this local static for correct initialization order
	if(dc) deprecatedMap[dc->name] = dc;
	return deprecatedMap;
}
map<string,DeprecatedCommand*>& getDeprecatedMap()
{	return updateDeprecatedMap(0);
}
DeprecatedCommand::DeprecatedCommand(string name) : name(name)
{	updateDeprecatedMap(this);
}


bool isReadable(string fname)
{	bool readable = false;
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname.c_str(), "r");
		if(fp)
		{	readable = true;
			fclose(fp);
		}
	}
	mpiWorld->bcast(readable);
	return readable;
}
