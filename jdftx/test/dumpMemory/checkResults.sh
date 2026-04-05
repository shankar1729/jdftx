#!/bin/bash

echo "4"  #number of checks

# Check 1: Reference SCF converges to expected energy
awk '/SCF: Cycle/ { E = $5 } END { print E, "-124.200 0.050 Reference SCF energy [Eh]" }' reference.out

# Check 2: dump-interval SCF converges to same energy as reference
# (tests that flushGpuCache does not corrupt computation)
Eref=$(awk '/SCF: Cycle/ { E = $5 } END { print E }' reference.out)
Edump=$(awk '/SCF: Cycle/ { E = $5 } END { print E }' withDumpInterval.out)
echo "$Edump $Eref 0.0001 dump-interval matches reference energy [Eh]"

# Check 3: dump-interval actually dumped intermediate files
nDumped=$(ls withDumpInterval.n* 2>/dev/null | wc -l)
if [ "$nDumped" -gt "0" ]; then dumped=1; else dumped=0; fi
echo "$dumped 1 0 dump-interval produced intermediate files"

# Check 4: .wfns file created by dump End State
if [ -f "withDumpInterval.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 dump End State created wfns file"
