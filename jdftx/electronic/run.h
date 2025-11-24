/*-------------------------------------------------------------------
Copyright 2025 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_RUN_H
#define JDFTX_ELECTRONIC_RUN_H

#include <core/string.h>

//Perform the actions of the main jdftx executable
//(except for the commandline/resource init and cleanup)
void run(string inputFilename, bool dryRun, bool printDefaults);

#endif //JDFTX_ELECTRONIC_RUN_H
