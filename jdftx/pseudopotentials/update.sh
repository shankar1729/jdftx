#!/bin/bash

DEST_DIR="$1"
SRC_DIR="$2"

#------ GBRV pseudoptentials ------
rm -rf $DEST_DIR/GBRV #Cleanup existing
tar -C $DEST_DIR -xpzf $SRC_DIR/GBRV.tgz
