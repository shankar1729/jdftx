#!/bin/bash

srcPath="$1"

echo '/** \page CommandIndex Index of commands'
cd "$srcPath/commands"
awk '$1=="/**" && $2=="\\page" {print "+ \\subpage ", $3}' *.dox | sort
echo '*/'
