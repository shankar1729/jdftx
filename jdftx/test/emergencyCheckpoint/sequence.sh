#!/bin/bash
#This test has a custom runner — do NOT list runs for the standard loop.
#Instead, we run jdftx ourselves, send SIGQUIT, and verify emergency dump.
export runs=""
export nProcs="1"

testSrcDir="$SRCDIR"
testRunDir="$(pwd)"
jdftxBuildDir="$4"

#Find the jdftx binary (passed as arg 4 to runTest.sh which sources us)
JDFTX="$jdftxBuildDir/jdftx$JDFTX_SUFFIX"
if [ ! -x "$JDFTX" ]; then
    #Fallback: try to find it from the positional args available in runTest.sh
    JDFTX="$(dirname "$(dirname "$testRunDir")")/jdftx$JDFTX_SUFFIX"
fi

echo "=== Emergency Checkpoint Test ==="
echo "Binary: $JDFTX"

#Clean previous results
rm -f emergency.out emergency.wfns emergency.force emergency.fillings emergency.eigenvals

#Run jdftx in background
$JDFTX -i "$testSrcDir/emergency.in" -d -o emergency.out &
JDFTX_PID=$!
echo "Started jdftx PID=$JDFTX_PID"

#Wait for SCF to begin (poll for first SCF line)
for i in $(seq 1 60); do
    if grep -q "SCF: Cycle" emergency.out 2>/dev/null; then
        echo "SCF started, sending SIGQUIT after brief delay..."
        sleep 1
        kill -QUIT $JDFTX_PID 2>/dev/null
        break
    fi
    sleep 1
done

#Wait for jdftx to exit cleanly
wait $JDFTX_PID 2>/dev/null
echo "jdftx exited with code $?"
