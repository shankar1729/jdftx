#!/bin/bash

echo "4"  #number of checks

# Check 1: WARNING present when no State dump configured
warnCount=$(grep -c "No 'dump.*State' command found" noStateDump.out 2>/dev/null)
echo "$warnCount 1 0 WARNING present when no dump State"

# Check 2: WARNING absent when State dump is configured
warnCount=$(grep -c "No 'dump.*State' command found" withStateDump.out 2>/dev/null)
echo "$warnCount 0 0 WARNING absent when dump End State"

# Check 3: .wfns file created by withStateDump
if [ -f "withStateDump.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 wfns file created with dump End State"

# Check 4: no .wfns file from noStateDump (no State in dump list, no SIGQUIT)
if [ -f "noStateDump.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 0 0 no wfns file without dump End State"
