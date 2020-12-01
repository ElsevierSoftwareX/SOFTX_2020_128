# -- doxygen -------

find_package(Doxygen)
find_package(LATEX COMPONENTS PDFLATEX MAKEINDEX DVIPDF BIBTEX)

if (DOXYGEN_FOUND)
  # set input and output files
  set(DOXYGEN_IN ./Doxyfile.in)
  set(DOXYGEN_OUT ./Doxyfile)
  
  # request to configure the file
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message("Doxygen build started")

  # copy auxiliary files
  file(
    COPY
    ./images
    ./style
    ./main.dox
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )
  
  # note the option ALL which allows to build the docs together with the application
  add_custom_target(doc_doxygen ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Doxygen"
    VERBATIM
    )

  # install doxygen
  install(
    DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
    )
  
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
