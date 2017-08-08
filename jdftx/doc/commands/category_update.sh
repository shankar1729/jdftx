#!/bin/bash

srcPath="$1"
cd "$srcPath/commands"

#Add input file format header:
echo '
/** \page Commands Input file documentation

Input file format
-----------------

+ The input file may contain commands in any order; commands will be
   automatically processed in an order that satisfies all dependencies

+ Each command is a single line, but lines can be continued using "\"

+ Whitespace is unimportant, except for separating words in each line

+ Any text following "#" on an input line is treated as a comment

+ "include @<file@>" can be used to include commands from @<file@>

+ Each instance of ${xxx} is replaced by environment variable "xxx"
   (Variable substitution occurs before command/include processing)

+ "set NAME VALUE" can be used to set an environment variable named NAME
   with value VALUE. This occurs before command/include processing,
   in the same pass as variable substitution. Therefore the order of
   "set xxx VALUE" and occurences of ${xxx} in the input file does matter

See \subpage CommandIndex for an alphabetical list of all available commands.

'

#Function to process for commands from a given executable
function processExecutable()
{
	#Heading for current executable
	if [ "$1" == "jdftx" ]; then
		echo "Commands allowed in jdftx input files"
	else
		echo "Additional commands in $1 input files"
	fi
	echo "-------------------------------------"
	
	awk '/^\/\/SectionInfo/ {print $2, $3, $4}' "$1.dox" | sort | awk '
	{	if(($1 != catPrev)    && ($1 != "NULL")) {cat=$1;    gsub("_"," ",cat);    printf("\n### %%%s\n", cat);}
		if(($2 != subCatPrev) && ($2 != "NULL")) {subCat=$2; gsub("_"," ",subCat); printf("\n<i>%%%s</i>\n", subCat);}
		printf("\\bigsep\\ref %s\n", $3);
		catPrev=$1; subCatPrev=$2;
	}
	'
	
	echo
}

processExecutable "jdftx"
processExecutable "phonon"
processExecutable "wannier"

#Close page
echo '*/'
