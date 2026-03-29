#!/bin/bash

echo "3"  #number of checks

# Check 1: SIGQUIT message present in output
sigquitMsg=$(grep -c "Received SIGQUIT" emergency.out 2>/dev/null)
echo "${sigquitMsg:-0} 1 0 SIGQUIT received and handled"

# Check 2: Emergency checkpoint message present
emergencyMsg=$(grep -c "Emergency checkpoint" emergency.out 2>/dev/null)
echo "${emergencyMsg:-0} 1 0 Emergency checkpoint message printed"

# Check 3: .wfns file was created despite no 'dump End State' in input
if [ -f "emergency.wfns" ]; then wfns=1; else wfns=0; fi
echo "$wfns 1 0 wfns file created by emergency dump"
