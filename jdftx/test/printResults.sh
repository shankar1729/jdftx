#!/bin/bash

srcDir="$1"
runDir="$2"

#Get list of test targets from CMakeLists.txt:
targets=( )
for target in $(grep '^add_jdftx_test' $srcDir/CMakeLists.txt); do
	target=${target/add_jdftx_test(/}
	target=${target/)/}
	targets+=( $target )
done

function printPaddedHeader()
{
	header="$1"
	netLen="78"
	headerLen="${#header}"
	leftPad=$(( ($netLen-$headerLen)/2 ))
	rightPad=$(( $netLen-$headerLen-$leftPad ))
	printf '=%.0s' $(eval echo "{1..$leftPad}")
	echo -n " ${header} "
	printf '=%.0s' $(eval echo "{1..$rightPad}")
	echo
}

#Print details of each test:
for target in "${targets[@]}"; do
	echo
	printPaddedHeader "Test ${target}"
	echo
	cat $runDir/$target/results
	cat $runDir/$target/summary
	echo
done

echo
printPaddedHeader "Summary of tests"
echo
for target in "${targets[@]}"; do
	printf "%30s: " ${target}
	cat $runDir/$target/summary
done
echo
