#!/bin/bash

echo "3"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.26856 0.0001 Vacuum energy [Eh]" }' vacuum.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.27986 0.0001 LinearPCM energy [Eh]" }' LinearPCM.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-17.28124 0.0001 CANDLE energy [Eh]" }' CANDLE.out
