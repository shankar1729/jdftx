#!/bin/bash

echo "4"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-125.95792 0.0001 TotalE Fe energy [Eh]" }' totalE.out
awk '/FillingsUpdate/ { M = $(NF-1) } END { print M, "+2.2983 0.010 TotalE Fe moment [muB]" }' totalE.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-125.95792 0.0001 SCF Fe energy [Eh]" }' SCF.out
awk '/FillingsUpdate/ { M = $(NF-1) } END { print M, "+2.2983 0.001 SCF Fe moment [muB]" }' SCF.out
