#!/bin/sh

sed -i 's|fftw3|fftw|g' ${OMICRONROOT}/${OMICRON_tag}/${OMICRONCONFIG}/build/omicron-*.*.*/src/CMakeLists.txt
