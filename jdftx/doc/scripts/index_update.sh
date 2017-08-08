#!/bin/bash

srcPath="$1"
cd "$srcPath/scripts"

scriptList=( ` find . -maxdepth 1 -type f -printf '%f\n' ` ) #get only files (ignore directories)

echo '/** \page Scripts Helper scripts

JDFTx includes several scripts to ease input generation and process / visualize output.
These scripts are located in the jdftx/scripts source directory.
Please add this directory to your PATH variable for convenience.
'

function getHeaderString()
{
	fname="$1"
	key="$2"
	awk "/$key/ {printf(\"%s\",\$2); for(i=3;i<=NF;i++) printf(\"_%s\",\$i);}" "$fname"
}

function sortedScriptList()
{
	for t in "${scriptList[@]}"; do
		echo "$(getHeaderString $t CATEGORY) $t $(getHeaderString $t SYNOPSIS)"
	done | sort
}

#---- Print links in main scripts page
sortedScriptList | awk '
	{	if($1 != catPrev) printf("\n+ \\subpage %s\n", $1);
		synopsis=$3; gsub("_", " ", synopsis);
		printf("  - \\ref %s : %s\n", $2, synopsis);
		catPrev=$1;
	}'

#---- Print links in sub-pages
sortedScriptList | awk '
	{	if($1 != catPrev)
		{	category=$1; gsub("_", " ", category);
			printf("*/\n\n/** \\page %s %s\n", $1, category);
		}
		synopsis=$3; gsub("_", " ", synopsis);
		printf("+ \\subpage %s : %s\n", $2, synopsis);
		catPrev=$1;
	}'
echo '*/'
echo

for t in "${scriptList[@]}"; do
	echo '/** \page ' $t $t
	echo "-----------------"
	./$t --help
	echo '*/'
	echo
done
