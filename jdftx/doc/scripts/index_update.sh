#!/bin/bash

srcPath="$1"

scriptList=( ` find "$srcPath/scripts" -maxdepth 1 -type f -printf '%f\n' ` ) #get only files (ignore directories)

echo '/** \page Scripts Scripts'
for t in "${scriptList[@]}"; do
	echo " + \subpage " $t
done
echo '*/'
echo

for t in "${scriptList[@]}"; do
	echo '/** \page ' $t $t
	echo "-----------------"
	$srcPath/scripts/$t --help
	echo '*/'
	echo
done
