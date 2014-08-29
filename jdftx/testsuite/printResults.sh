#!/bin/bash

cd $1
shopt -s nullglob

for file in */results; do
	echo "Test ${file/\/results/}:"
	cat $file
	cat ${file/results/summary}
	echo
done

echo "Summary of tests:"
for file in */summary; do
	printf "%20s: " ${file/\/summary/}
	cat $file
done
echo