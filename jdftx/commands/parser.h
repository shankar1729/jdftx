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

#ifndef JDFTX_PARSER_H
#define JDFTX_PARSER_H

#include <electronic/common.h>

//! Parse the input file, and initialize everything
//! Read from stdin if filename is empty
void parse(const char* filename, Everything &everything, bool printDefaults=false);

//! Produce an input file template
void printDefaultTemplate(Everything &e);

//! Write an HTML manual to stdout
void writeCommandManual(Everything& e);

#endif // JDFTX_PARSER_H
