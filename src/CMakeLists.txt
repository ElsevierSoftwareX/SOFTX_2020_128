# -- dependencies -----------

# FrameL
pkg_check_modules(FRAMEL REQUIRED "libframe")

# HDF5
find_package(HDF5)

# ROOT
find_package(ROOT)
message(STATUS "Found ROOT libraries in ${ROOT_LIBRARY_DIR}")
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
include(${ROOT_USE_FILE})

# GWOLLUM
#find_package(GWOLLUM 2.3.11 REQUIRED)
pkg_check_modules(GWOLLUM REQUIRED IMPORTED_TARGET "gwollum")
message(STATUS "GWOLLUM headers: ${GWOLLUM_INCLUDE_DIRS}")
message(STATUS "GWOLLUM libraries: ${GWOLLUM_LIBRARY_DIRS}")
#message(STATUS "GWOLLUM libraries: ${GWOLLUM_LIBRARIES}")
message(STATUS "HDF5 libraries: ${HDF5_INCLUDE_DIRS}")

# add link paths
link_directories(${FRAMEL_LIBRARY_DIRS} ${GWOLLUM_LIBRARY_DIRS})

# versioning
configure_file(Oconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/Oconfig.h @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


# -- libOmicron -------------

set(
  OMICRON_HEADERS
  Omap.h
  Oqplane.h
  Otile.h
  Oinject.h
  Oomicron.h
  )

# build ROOT dictionary
include_directories(
  ${GWOLLUM_INCLUDE_DIRS}
  ${FRAMEL_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  )
root_generate_dictionary(
  OmicronDict
  LINKDEF LinkDef.h
  MODULE Omicron
  OPTIONS ${OMICRON_HEADERS} -I${CMAKE_CURRENT_SOURCE_DIR}
  )

# compile library
# NOTE: we call this `libOmicron` here to not clash with the `omicron`
#       executable target below, cmake can get easily confused on
#       case-insensitive file systems
add_library(
  libOmicron
  SHARED
  Ohtml.cc
  Oinject.cc
  Omap.cc
  Oomicron.cc
  OmicronDict.cxx
  Oparameters.cc
  Oqplane.cc
  Otile.cc
  Outils.cc
  )
target_include_directories(
  libOmicron
  PUBLIC
  ${GWOLLUM_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(
  libOmicron
  ${ROOT_Core_LIBRARY}
  ${ROOT_Hist_LIBRARY}
  ${ROOT_MathCore_LIBRARY}
  ${ROOT_RIO_LIBRARY}
  CUtils
  Inject
  RootUtils
  Segments
  Time
  Triggers
  )
set_target_properties(
  libOmicron PROPERTIES
  OUTPUT_NAME Omicron
  )

# -- libOmicronUtils --------

add_library(
  OmicronUtils
  SHARED
  OmicronUtils.cc
  )
target_link_libraries(
  OmicronUtils
  CUtils
  Streams
  )

# -- executables ------------

add_executable(
  omicron
  omicron.cc
  )
target_link_libraries(
  omicron
  ${ROOT_Core_LIBRARY}
  CUtils
  Segments
  libOmicron
  )

add_executable(
  omicron-listfile
  omicron-listfile.cc
  )
target_link_libraries(
  omicron-listfile
  OmicronUtils
  )

add_executable(
  omicron-scanfile
  omicron-scanfile.cc
  )
target_link_libraries(
  omicron-scanfile
  ${ROOT_Core_LIBRARY}
  ${ROOT_RIO_LIBRARY}
  CUtils
  OmicronUtils
  )

add_executable(
  omicron-plot
  omicron-plot.cc
  )
target_link_libraries(
  omicron-plot
  ${ROOT_Core_LIBRARY}
  ${ROOT_Graf_LIBRARY}
  CUtils
  RootUtils
  Segments
  Triggers
  OmicronUtils
  )

add_executable(
  omicron-print
  omicron-print.cc
  )
target_link_libraries(
  omicron-print
  ${ROOT_Core_LIBRARY}
  CUtils
  Segments
  Triggers
  OmicronUtils
  )

# -- installation -----------


# install library
install(
  TARGETS
  libOmicron
  OmicronUtils
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# install ROOT PCM
install(
  FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libOmicron_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/libOmicron.rootmap
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

# install executable(s)
install(
  TARGETS
  omicron
  omicron-listfile
  omicron-scanfile
  omicron-plot
  omicron-print
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# -- pkgconfig --------------

# install pkgconfig file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_LOWER}.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}.pc"
  @ONLY
  )
install(
  FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}.pc
  DESTINATION
  ${CMAKE_INSTALL_LIBDIR}/pkgconfig/
  )
