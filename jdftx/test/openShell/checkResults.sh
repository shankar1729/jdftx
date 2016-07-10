#!/bin/bash

echo "4"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-0.50044 0.0001 H energy [Eh]" }' Hatom.out
awk '/magnetic-moments H/ { M = $NF } END { print M, "+0.997 0.010 H magnetic moment [muB]" }' Hatom.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-32.02300 0.0001 O2 energy [Eh]" }' O2.out
awk '/magnetic-moments O/ { M = $NF } END { print M, "+0.995 0.010 O magnetic moment in O2 [muB]" }' O2.out
