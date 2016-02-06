#!/bin/bash
major=$(grep "set(CPACK_PACKAGE_VERSION_MAJOR" jdftx/CMakeLists.txt | cut --delimiter="\"" -f 2)
minor=$(grep "set(CPACK_PACKAGE_VERSION_MINOR" jdftx/CMakeLists.txt | cut --delimiter="\"" -f 2)
patch=$(grep "set(CPACK_PACKAGE_VERSION_PATCH" jdftx/CMakeLists.txt | cut --delimiter="\"" -f 2)
filename="jdftx.$major.$minor.$patch"

tar -cpzf "../$filename.tar.gz" jdftx --exclude='.*' --exclude='*~' --exclude='*.tgz'

