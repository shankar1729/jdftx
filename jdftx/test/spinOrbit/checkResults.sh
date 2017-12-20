#!/bin/bash

echo "2"  #number of checks

awk '/IonicMinimize: Iter/ { E = $5 } END { print E, "-88.9619 0.0001 Pt bulk energy [Eh]" }' Pt.out
awk '/FillingsUpdate/ { mu = $3 } END { print mu, "0.8333 0.0001 Pt bulk mu [Eh]" }' Pt.out
