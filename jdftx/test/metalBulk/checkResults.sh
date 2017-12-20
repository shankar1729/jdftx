#!/bin/bash

echo "4"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-124.3922 0.0001 TotalE Fe energy [Eh]" }' totalE.out
awk '/FillingsUpdate/ { M = $(NF-1) } END { print M, "+2.272 0.001 TotalE Fe moment [muB]" }' totalE.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-124.3922 0.0001 SCF Fe energy [Eh]" }' SCF.out
awk '/FillingsUpdate/ { M = $(NF-1) } END { print M, "+2.272 0.001 SCF Fe moment [muB]" }' SCF.out
