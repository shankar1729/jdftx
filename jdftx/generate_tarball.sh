version_major=$(grep "set(CPACK_PACKAGE_VERSION_MAJOR" CMakeLists.txt | cut --delimiter="\"" -f 2)
version_minor=$(grep "set(CPACK_PACKAGE_VERSION_MINOR" CMakeLists.txt | cut --delimiter="\"" -f 2)
version_patch=$(grep "set(CPACK_PACKAGE_VERSION_PATCH" CMakeLists.txt | cut --delimiter="\"" -f 2)
filename="jdftx.$version_major.$version_minor.$version_patch"

tar -c --gzip -f "../$filename.tar.gz" * --exclude='.*' --exclude='*~'
