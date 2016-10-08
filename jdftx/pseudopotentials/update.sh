#!/bin/bash

DEST_DIR="$1"
SRC_DIR="$2"
URL_PREFIX="http://sourceforge.net/p/jdftx/git/ci/HEAD/tree/jdftx/pseudopotentials/"

function processSet()
{
	setName="$1"
	srcFile="$SRC_DIR/${setName}.tgz"
	if [ ! -f $srcFile ]; then
		wget -O "$srcFile" "${URL_PREFIX}${setName}.tgz?format=raw"
	fi
	rm -rf "$DEST_DIR/$setName" #Cleanup existing
	tar -C "$DEST_DIR" -xpzf "$srcFile"

}

processSet GBRV
processSet SG15
processSet SG15-pulay

