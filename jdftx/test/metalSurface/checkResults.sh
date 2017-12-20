#!/bin/bash

echo "4"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-266.8440 0.0001 Neutral Pt(111) energy [Eh]" }' neutral.out
awk '/FillingsUpdate/ { mu = $3 } END { print mu, "-0.196 0.001 Neutral Pt(111) mu [Eh]" }' neutral.out
awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-259.1668 0.0001 Fixed-mu Pt(111) energy [Eh]" }' charged.out
awk '/FillingsUpdate/ { Q = $NF - 48 } END { print Q, "0.211 0.001 Fixed-mu Pt(111) charge [e]" }' charged.out
