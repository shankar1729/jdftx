#!/bin/bash

echo "3"  #number of checks

awk '/Single-point solvation/ { print $NF, "-0.07372 0.0001 Single-point solvation [Eh]" }' LinearPCM.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-24.47307 0.0001 LinearPCM energy [Eh]" }' LinearPCM.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-24.51048 0.0001 CANDLE energy [Eh]" }' CANDLE.out
