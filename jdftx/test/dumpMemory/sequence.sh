#!/bin/bash
#Test that dump-interval does not affect SCF convergence (memory leak fix)
#and that out-of-core Pulay history produces same results as baseline.
export runs="reference withDumpInterval"
export nProcs="1"
