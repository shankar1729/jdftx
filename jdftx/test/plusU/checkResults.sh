#!/bin/bash

echo "2"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-375.25826 0.0001 NiO AFM energy [Eh]" }' NiO.out
awk '/FillingsUpdate/ { M = $(NF-3) } END { print M, "3.49624 0.001 NiO AFM abs moment [muB]" }' NiO.out
