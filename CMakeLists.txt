#
# CMake packaging for Omicron
# Copyright 2019  Duncan Macleod <duncan.macleod@ligo.org>
#

# -- package info -----------

cmake_minimum_required(VERSION 3.9.0 FATAL_ERROR)
project(
    Omicron
    LANGUAGES CXX
    VERSION 2.3.11
    DESCRIPTION "An algorithm to detect and characterize transient events in gravitational-wave detectors"
)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

# -- extra requirements -----

include(GNUInstallDirs)
find_package(PkgConfig REQUIRED)

# -- set install paths ------

if (WIN32)
    set(RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
    set(LIBRARY_DESTINATION ${CMAKE_INSTALL_LIBDIR})
else (WIN32)
    set(RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
    set(LIBRARY_DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(WIN32)

# -- build components -------

add_subdirectory(src)

# -- build tarball ----------
#
# to build a source tarball:
#
# mkdir dist
# pushd dist
# cmake ..
# cmake --build . --target package_source
#
# this should produce dist/omicron-X.Y.Z.tar.xz
#

set(CPACK_PACKAGE_MAJOR ${${PROJECT_NAME}_MAJOR_VERSION})
set(CPACK_PACKAGE_MINOR ${${PROJECT_NAME}_MINOR_VERSION})
set(CPACK_PACKAGE_PATCH ${${PROJECT_NAME}_PATCH_VERSION})
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_MAJOR}.${CPACK_PACKAGE_MINOR}.${CPACK_PACKAGE_PATCH}")

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

set(CPACK_SOURCE_GENERATOR TXZ)
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROJECT_NAME_LOWER}-${${PROJECT_NAME}_VERSION})
set(CPACK_SOURCE_IGNORE_FILES
    "/.*~$/"
    ".*~$"
    "\\\\.deps/"
    "\\\\.git/"
    "\\\\.svn/"
    "CMakeFiles/"
    "_CPack_Packages/"
    "build/"
    "doc/"
    "html/"
    "paper/"
    "\\\\.cmake"
    "\\\\.tar\\\\.xz"
    "CMakeCache.txt"
    "Makefile"
)
include(CPack)
