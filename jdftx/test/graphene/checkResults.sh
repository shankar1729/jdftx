#!/bin/bash

echo "4"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-11.42190 0.0001 Graphene energy [Eh]" }' totalE.out
awk '/FillingsUpdate/ { mu = $3 } END { print mu, "-0.15056 0.0001 Graphene mu [Eh]" }' totalE.out

which octave > /dev/null
if [ $? -eq 0 ]; then
	${SRCDIR}/../../scripts/binaryToText bandstruct.eigenvals > bandstruct.eigenvals.txt
	awk 'NF==8 { print $4, "-0.15056 0.0001 Dirac point band 1 [Eh]" }' bandstruct.eigenvals.txt
	awk 'NF==8 { print $5, "-0.15056 0.0001 Dirac point band 2 [Eh]" }' bandstruct.eigenvals.txt
else
	echo "0 0 1 octave not found"
	echo "0 0 1 octave not found"
fi
