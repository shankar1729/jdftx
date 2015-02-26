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
{	static map<string,Command*> commandMap; //make this local static for correct initializtaion order
	if(command) commandMap[command->name] = command;
	return commandMap;
}
map<string,Command*>& getCommandMap()
{	return updateCommandMap(0);
}

Command::Command(string name, string section)
: name(name), section(section), allowMultiple(false), hasDefault(false)
{	updateCommandMap(this);
}

void Command::require(string name) { requires.insert(name); }
void Command::forbid(string name) { forbids.insert(name); }

bool isReadable(string fname)
{	bool readable = false;
	if(mpiUtil->isHead())
	{	FILE* fp = fopen(fname.c_str(), "r");
		if(fp)
		{	readable = true;
			fclose(fp);
		}
	}
	mpiUtil->bcast(readable);
	return readable;
}
