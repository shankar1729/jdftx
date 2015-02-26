#!/bin/bash

srcPath="$1"
outFile="$2"

echo '/** \page CommandIndex Index of commands' > $outFile
awk '$1=="/**" && $2=="\\page" {print "+ \\subpage ", $3}' $srcPath/*/manual.dox | sort >> $outFile
echo '*/' >> $outFile
