# add a target to generate API documentation with Doxygen
find_package (Doxygen)

if (${DOXYGEN_FOUND})

  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doxy-codeare.conf.in
    ${CMAKE_CURRENT_BINARY_DIR}/doxy-codeare.conf @ONLY) 

  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxy-codeare.conf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/codeare)

endif()
