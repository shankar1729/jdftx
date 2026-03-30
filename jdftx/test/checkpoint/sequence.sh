#!/bin/bash
#Combined checkpoint test: dump State, OMP_NUM_THREADS, and SIGQUIT emergency dump.
#Uses custom runner (runs="") since SIGQUIT test needs background process control.
export runs=""
export nProcs="1"

export OMP_NUM_THREADS=3

testSrcDir="$SRCDIR"
JDFTX="$4/jdftx$JDFTX_SUFFIX"
if [ ! -x "$JDFTX" ]; then
    JDFTX="$(dirname "$(dirname "$(pwd)")")/jdftx$JDFTX_SUFFIX"
fi

#Run 1: Standard run with dump End State (also verifies OMP_NUM_THREADS)
echo "=== Run: withStateDump (OMP_NUM_THREADS=$OMP_NUM_THREADS) ==="
$JDFTX -i "$testSrcDir/withStateDump.in" -d -o withStateDump.out

#Run 2: Emergency checkpoint via SIGQUIT
echo "=== Run: emergency SIGQUIT test ==="
rm -f emergency.wfns
$JDFTX -i "$testSrcDir/emergency.in" -d -o emergency.out &
JDFTX_PID=$!
echo "Started jdftx PID=$JDFTX_PID"

for i in $(seq 1 60); do
    if grep -q "SCF: Cycle" emergency.out 2>/dev/null; then
        echo "SCF started, sending SIGQUIT..."
        sleep 1
        kill -QUIT $JDFTX_PID 2>/dev/null
        break
    fi
    sleep 1
done

wait $JDFTX_PID 2>/dev/null
echo "jdftx exited with code $?"
