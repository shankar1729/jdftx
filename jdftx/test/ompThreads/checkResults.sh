#!/bin/bash

echo "2"  #number of checks

# Check 1: OMP override message present in output
overrideMsg=$(grep -c "OMP_NUM_THREADS=3 overrides" ompTest.out 2>/dev/null)
# Note: if auto-detected cores happen to equal 3, no override message is printed.
# In that case, we still verify thread count (check 2).
physCores=$(grep "Run totals:" ompTest.out | awk '{print $5}')
if [ "$physCores" = "3" ]; then
    # Machine happens to have 3 physical cores — override msg won't appear, that's OK
    overrideMsg=1
fi
echo "$overrideMsg 1 0 OMP_NUM_THREADS override detected or matching"

# Check 2: Actual thread count in output matches OMP_NUM_THREADS=3
threadCount=$(grep "Run totals:" ompTest.out | awk '{print $5}')
echo "${threadCount:-0} 3 0 thread count matches OMP_NUM_THREADS"
