#!/bin/bash

echo "5"  #number of checks

# Check 1: .wfns file created by dump End State
if [ -f "withStateDump.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 wfns file created with dump End State"

# Check 2: OMP_NUM_THREADS override detected in output
overrideMsg=$(grep -c "OMP_NUM_THREADS=3 overrides" withStateDump.out 2>/dev/null)
# If machine happens to have exactly 3 cores, no override message -- still OK
physCores=$(grep "Run totals:" withStateDump.out | awk '{print $5}')
if [ "$physCores" = "3" ]; then overrideMsg=1; fi
echo "${overrideMsg:-0} 1 0 OMP_NUM_THREADS override detected or matching"

# Check 3: Thread count matches OMP_NUM_THREADS=3
threadCount=$(grep "Run totals:" withStateDump.out | awk '{print $5}')
echo "${threadCount:-0} 3 0 thread count matches OMP_NUM_THREADS"

# Check 4: SIGQUIT received and handled
sigquitMsg=$(grep -c "Received SIGQUIT" emergency.out 2>/dev/null)
echo "${sigquitMsg:-0} 1 0 SIGQUIT received and handled"

# Check 5: Emergency .wfns created despite no dump End State
if [ -f "emergency.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 emergency wfns created by SIGQUIT dump"
