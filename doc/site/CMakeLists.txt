# add a target to generate API documentation with Doxygen
find_package (XSLTPROC)

if (${XSLTPROC_FOUND})

  list (APPEND XMLPAGES xml/contact.xml xml/developerguide.xml
    xml/download.xml xml/gettingstarted.xml xml/index.xml
    xml/install.xml xml/userguide.xml)  

  set (page_xsl ${CMAKE_CURRENT_SOURCE_DIR}/xsl/codeare.org.xsl)
  
  set(htmlfiles)
  foreach(page_xml ${XMLPAGES})
    string(REPLACE ".xml" ".html" page_html ${page_xml})
    string(REPLACE "xml/" "" page_html ${page_html})
    add_custom_command (OUTPUT ${page_html}
      COMMAND "${XSLTPROC_EXECUTABLE}"
      ARGS ${page_xsl}
        ${CMAKE_CURRENT_SOURCE_DIR}/${page_xml} >
        ${CMAKE_CURRENT_BINARY_DIR}/${page_html} 
      DEPENDS ${page_xml} ${page_xsl}
      COMMENT "Regenerating ${page_html} with XSLTPROC" VERBATIM
      )
    string (REPLACE ".html" "" page_target ${page_html})
    add_custom_target (${page_html} DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${page_html})
    list (APPEND htmlfiles ${CMAKE_CURRENT_BINARY_DIR}/${page_html})
    list (APPEND htmltargets ${page_target})
  endforeach()

  add_custom_target (site DEPENDS ${htmlfiles})

  list (APPEND jsfiles nav.js sh_cpp.js sh_main.js sh_sh.js sh_xml.js)
  list (APPEND cssfiles screen.css sh_style.css)
  install (FILES ${cssfiles} ${htmlfiles} ${jsfiles}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/codeare) 

endif()
