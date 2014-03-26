
if (NOT XSLTPROC_FOUND)
  IF (WIN32)
    SET (xsltproc_n "xsltproc.exe")
  ELSE()
    SET (xsltproc_n "xsltproc")
  ENDIF()
  
  FIND_PROGRAM(XSLTPROC_EXECUTABLE ${xsltproc_n} PATHS "${EXTERNALS_DIRECTORY}/xsltproc/${osdir}")
  
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(XSLTPROC DEFAULT_MSG
    XSLTPROC_EXECUTABLE
    )
  MARK_AS_ADVANCED(XSLTPROC_EXECUTABLE)
ENDIF()
