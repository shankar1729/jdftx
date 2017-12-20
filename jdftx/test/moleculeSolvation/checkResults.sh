#!/bin/bash

echo "3"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.2681 0.0001 Vacuum energy [Eh]" }' vacuum.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.2793 0.0001 LinearPCM energy [Eh]" }' LinearPCM.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.2807 0.0001 CANDLE energy [Eh]" }' CANDLE.out
