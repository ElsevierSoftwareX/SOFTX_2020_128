#!/bin/bash

# get cmake parameters
PROJECT_NAME=`grep -m1 -A1 "project(" ./CMakeLists.txt | grep -v "project("`
VERSION=`grep -m1 " VERSION " ./CMakeLists.txt | awk '{print $2}' | sed 's|\.|\\\.|g'`
DESCRIPTION=`grep -m1 " DESCRIPTION " ./CMakeLists.txt | cut -d '"' -f2`

# generate doxygen configuration
sed -e "s|@PROJECT_NAME@|$PROJECT_NAME|g" \
    -e "s|@PROJECT_VERSION@|${VERSION}|g" \
    -e "s|@PROJECT_DESCRIPTION@|${DESCRIPTION}|g" \
    -e "s|@CMAKE_CURRENT_BINARY_DIR@|\./doc|g" \
    -e "s|@CMAKE_SOURCE_DIR@|\.|g" \
    ./doc/Doxyfile.in > ./doc/Doxyfile

# generate html doc
doxygen ./doc/Doxyfile

