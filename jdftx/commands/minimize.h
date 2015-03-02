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

#ifndef JDFTX_COMMAND_MINIMIZE_H
#define JDFTX_COMMAND_MINIMIZE_H

/** @file minimize.h
@brief Provides base class for defining the minimize commands
*/

#include <commands/command.h>
#include <electronic/Everything.h>

//! @brief Abstract base class for all the minimize commands
struct CommandMinimize : public Command
{	CommandMinimize(string systemName, string section="jdftx"); //!< provide a command called systemName-minimize in specified section
	void process(ParamList& pl, Everything& e);
	void printStatus(Everything& e, int iRep);
protected:
	//! @brief Derived class should specify where the parameters are stored
	//! @param e Reference to Everything
	//! @return Return reference to the relevant MinimizeParams that this command should operate on.
	virtual MinimizeParams& target(Everything& e)=0;
};

#endif // JDFTX_COMMAND_MINIMIZE_H
