#!/bin/bash

DEST_DIR="$1"
SRC_DIR="$2"

#------ GBRV pseudoptentials ------
GBRV_PATH="$DEST_DIR/GBRV"

echo "Downloading GBRV pseudopotentials from http://www.physics.rutgers.edu/gbrv to $GBRV_PATH
     K.F. Garrity, J.W. Bennett, K.M. Rabe and D. Vanderbilt, Comput. Mater. Sci. 81, 446 (2014)
These pseudopotentials are distributed under the GPLv3 license. Please see the
above website for copyright information as well as details on the pseudopotentials."

mkdir -p $GBRV_PATH
function downloadGBRV()
{	wget --timestamping --max-redirect=0 --no-verbose --directory-prefix=$GBRV_PATH http://www.physics.rutgers.edu/gbrv/$1
	tar -C $GBRV_PATH -xpzf $GBRV_PATH/$1
}
downloadGBRV all_pbe_uspp.tar.gz
downloadGBRV all_lda_uspp.tar.gz
