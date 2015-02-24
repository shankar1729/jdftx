#!/bin/bash
for w in 55 128 256 384; do
	dpi=$(echo $w | awk '{print $w*90.0/128}')
	inkscape -e jdftx-$w.png -d $dpi -C jdftx.svg
done