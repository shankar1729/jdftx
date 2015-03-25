#!/bin/bash

srcPath="$1"

echo '/** \page CommandIndex Index of commands'
awk '$1=="/**" && $2=="\\page" {print "+ \\subpage ", $3}' $srcPath/commands/*.dox | sort
echo '*/'
