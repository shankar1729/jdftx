#!/bin/bash

srcPath="$1"
outPath="$2"

echo '/** \page Scripts Scripts' > $outPath/Scripts.dox

scriptList=( `ls "$srcPath/scripts" `)

 for t in "${scriptList[@]}"
do
if [ "$t" = "ase" ]; then
  echo "Skipping ase directory"
else
   echo " + \subpage " $t >> $outPath/Scripts.dox
   echo '/** \page ' $t $t > $outPath/$t.dox    
   echo "Details and usage:" >> $outPath/$t.dox
     if [ "$t" = "dryRunToPDB" ]; then
       $srcPath/scripts/$t >> $outPath/$t.dox
     else
       $srcPath/scripts/$t --help >> $outPath/$t.dox
  fi
   echo '*/' >> $outPath/$t.dox
fi
done

echo '*/' >> $outPath/Scripts.dox
