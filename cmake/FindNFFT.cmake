# - Find NFFT library
# Find the native NFFT includes and library
# This module defines
#  NFFT_INCLUDE_DIR, where to find hdf5.h, etc.
#  NFFT_LIBRARIES, libraries to link against to use NFFT.
#  NFFT_FOUND, If false, do not try to use NFFT.
# also defined, but not for general use are
#  NFFT_LIBRARY, where to find the NFFT library.

FIND_PATH(NFFT_INCLUDE_DIR nfft3.h PATH_SUFFIXES nfft)

SET(NFFT_NAMES ${NFFT_NAMES} nfft nfft3)
FIND_LIBRARY(NFFT_LIBRARY NAMES ${NFFT_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set NFFT_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NFFT DEFAULT_MSG NFFT_LIBRARY NFFT_INCLUDE_DIR)

IF(NFFT_FOUND)
  SET( NFFT_LIBRARIES ${NFFT_LIBRARY} )
ENDIF(NFFT_FOUND)
