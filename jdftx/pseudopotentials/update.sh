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
processSet GBRV-2017
processSet SG15
processSet SG15-pulay


#Add links to latest GBRV pseudopotentials:
cd $DEST_DIR/GBRV
prefixes=$( for file in *_*_*.uspp; do
		basename $file | awk '{split($1,tokens,"_v"); print tokens[1]}'
	done | sort | uniq )
for prefix in $prefixes; do
	for version in "1" "1.01" "1.2" "1.4" "1.5"; do
		fname="${prefix}_v${version}.uspp"
		if [ -f $fname ]; then
			latestVersion="$fname"
		fi
	done
	ln -sf $latestVersion ${prefix}.uspp
done

#Add links to latest SG15 pseudopotentials:
cd $DEST_DIR/SG15
prefixes=$( for file in *_*_*-*.upf; do
		basename $file | awk '{split($1,tokens,"-"); print tokens[1]}'
	done | sort | uniq )
for prefix in $prefixes; do
	for suffix in "upf" "pulay"; do
		for version in "1.0" "1.1"; do
			fname="${prefix}-${version}.${suffix}"
			if [ -f $fname ]; then
				latestVersion="$fname"
			fi
		done
		ln -sf $latestVersion ${prefix}.${suffix}
	done
done
