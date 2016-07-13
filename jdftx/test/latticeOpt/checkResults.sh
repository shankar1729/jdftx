#!/bin/bash

echo "3"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-7.93641 0.0001 Si energy [Eh]" }' latticeOpt.out
awk 'NR==2 {x0=$3} NR==3 { print $3-x0, "0.25 0.001 Si fractional coordinate" }' latticeOpt.ionpos
awk 'NR==2 { print sqrt($1*$1+$2*$2+$3*$3), "7.3216 0.001 Si latvec length [a0]" }' latticeOpt.lattice
