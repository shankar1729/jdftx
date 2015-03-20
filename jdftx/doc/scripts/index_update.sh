#!/bin/bash

srcPath="$1"
outFile="$2"

echo '/** \page Scripts Scripts' > $outFile

scriptList=( `ls "$srcPath/scripts" `)

 for t in "${scriptList[@]}"
do
if [ "$t" = "ase" ]; then
  echo "Skipping ase directory"
else
  if [ "$t" = "dryRunToPDB" ]; then
    echo $t >> $outFile
    $srcPath/scripts/$t >> $outFile
  else
    echo $t >> $outFile
    $srcPath/scripts/$t --help >> $outFile
  fi
fi
done

echo '*/' >> $outFile



