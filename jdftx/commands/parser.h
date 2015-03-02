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

/** @file parser.h
@brief Functions for parsing JDFTx input files
*/

#include <electronic/common.h>
#include <vector>

//! @brief Read input from file into an array of command/argument-list pairs.
//! This routine handles comments, line continuation as well as environment substitution.
//! In MPI mode, the input will be read only on head and distributed to all the processes.
//! @param filename Name of file to read, or empty for reading from stdin
//! @return Pairs of commands and their argument lists
std::vector< std::pair<string,string> > readInputFile(string filename);

//! @brief Parse input and initialize everything
//! @param input List o commands and arguments as obtained by readInputFile()
//! @param everything Reference to everything to initialize
//! @param printDefaults If true, print status of commands that were invoked by default in addition to those invoked manually
void parse(std::vector< std::pair<string,string> > input, Everything &everything, bool printDefaults=false);

//! @brief Output command documentation in the form of an input file template.
//! This will output the comments for all the commands, along with the default command line, if any.
//! @param e A reference to Everything
void printDefaultTemplate(Everything &e);

//! @brief Output command docuemntation in Doxygen format which gets converted to the HTML and LaTeX manuals
//! @param e A reference to Everything
//! @param section Which executable to generate documentation for i.e. only include commands whose section matches this
void writeCommandManual(Everything& e, string section);

#endif // JDFTX_PARSER_H
