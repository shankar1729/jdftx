#!/bin/bash

if [ -z "$1" ]; then echo "First command-line argument should be absolute location of xc_funcs.h from LibXC"; exit 1; fi

processedData="$(awk '/XC_/ {
	xcName = tolower(substr($2,4));
	gsub(/_/, "-", xcName);
	printf("%-30s, \"%s\",\n", $2, xcName); }' $1)"

echo "$processedData" | grep "XC_LDA_X "
echo "$processedData" | grep "XC_LDA_X_";
echo
echo "$processedData" | grep "XC_GGA_X_"
echo
echo "$processedData" | grep "XC_MGGA_X_"
echo
echo
echo "$processedData" | grep "XC_LDA_C_";
echo
echo "$processedData" | grep "XC_GGA_C_"
echo
echo "$processedData" | grep "XC_MGGA_C_"
echo
echo
echo "$processedData" | grep "XC_LDA_XC_";
echo
echo "$processedData" | grep "XC_GGA_XC_"
echo
echo "$processedData" | grep "XC_MGGA_XC_"
echo
echo "$processedData" | grep "XC_HYB_GGA_XC_"
echo
echo "$processedData" | grep "XC_HYB_MGGA_XC_"
echo
echo
echo "$processedData" | grep "XC_LDA_K_"
echo
echo "$processedData" | grep "XC_GGA_K_"
echo
echo
echo

