#!/bin/bash

srcPath="$1"
cd "$srcPath/scripts"

scriptList=( ` find . -maxdepth 1 -type f -printf '%f\n' ` ) #get only files (ignore directories)

echo '/** \page Scripts Scripts

JDFTx includes several scripts to ease input generation and process / visualize output.
These scripts are located in the jdftx/scripts source directory.
Please add this directory to your PATH variable for convenience.

'
for t in "${scriptList[@]}"; do
	echo " + \subpage " $t
done
echo '*/'
echo

for t in "${scriptList[@]}"; do
	echo '/** \page ' $t $t
	echo "-----------------"
	./$t --help
	echo '*/'
	echo
done
