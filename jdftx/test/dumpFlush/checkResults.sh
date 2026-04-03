#!/bin/bash

echo "4"  #number of checks

# Check 1: Reference SCF converges
awk '/SCF: Cycle/ { E = $5 } END { print E, "-124.200 0.050 Reference SCF energy [Eh]" }' reference.out

# Check 2: dump-interval produces same energy (cache flush does not corrupt)
Eref=$(awk '/SCF: Cycle/ { E = $5 } END { print E }' reference.out)
Edump=$(awk '/SCF: Cycle/ { E = $5 } END { print E }' withDumpInterval.out)
echo "$Edump $Eref 0.0001 dump-interval matches reference energy [Eh]"

# Check 3: Intermediate density files were actually dumped
nDumped=$(ls withDumpInterval.n* 2>/dev/null | wc -l)
if [ "$nDumped" -gt "0" ]; then dumped=1; else dumped=0; fi
echo "$dumped 1 0 dump-interval produced intermediate density files"

# Check 4: Final state was dumped
if [ -f "withDumpInterval.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 dump End State created wfns file"
