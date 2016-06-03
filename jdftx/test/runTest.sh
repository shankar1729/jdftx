#!/bin/bash

testName="$1"
testsuiteSrcDir="$2"
testsuiteRunDir="$3"
jdftxBuildDir="$4"

testSrcDir="$testsuiteSrcDir/$testName"
testRunDir="$testsuiteRunDir/$testName"

mkdir -p $testRunDir
cd $testRunDir

#Run JDFTx on all the runs that belong to this test (don't rerun tests which have succeeded)
source $testSrcDir/sequence.sh
for run in $runs; do
	if [[ ! ( ( -f $run.out ) && ( "$(awk '/End date and time:/ {endLine=NR+1} NR==endLine {print}' $run.out)" == "Done!" ) ) ]]; then
		export SRCDIR="$testSrcDir"
		$JDFTX_LAUNCH $jdftxBuildDir/jdftx -i $testSrcDir/$run.in -d -o $run.out
		if [ "$?" -ne "0" ]; then
			echo "" > results
			echo "FAILED: error running $run" > summary
			exit 1
		fi
	fi
done

#Parse results:
$testSrcDir/checkResults.sh | awk '
	NR==1 { nChecks = $1; iCheck = 0; nFail = 0; 
		printf("%30s  %19s %19s Status\n", "Check name", "Obtained value", "Expected value");
	}
	NR>1 && NF>4 {
		iCheck++;
		xObtained = $1;
		xExpected = $2;
		xTolerance = $3;
		checkName = $4; for(k=5; k<=NF; k++) checkName = checkName " " $k;
		fail = (xObtained < xExpected-xTolerance) || (xObtained > xExpected+xTolerance);
		if(fail) nFail++;
		printf("%30s: %19.12le %19.12le [%s]\n", checkName, xObtained, xExpected, fail ? "FAILED" : "Passed");
	}
	END {
		if(iCheck != nChecks)
		{	print "FAILED: mismatch in number of checks (most likely a parse error)" > "summary";
			exit 1;
		}
		else if(nFail > 0)
		{	printf("FAILED: failed %d of %d checks.\n", nFail, nChecks) > "summary";
			exit 1;
		}
		else
		{	printf("Passed: passed %d checks.\n", nChecks) > "summary";
			exit 0;
		}
	}
' > results
