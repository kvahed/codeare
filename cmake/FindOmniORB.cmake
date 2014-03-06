#
# Find the omniORB libraries and include dir
#
 
# The following variables are set:
# OMNIORB4_FOUND        - When false, don't try to use omniORB
# OMNIORB4_INCLUDE_DIR  - Directories to include to use omniORB
# OMNIORB4_LIBRARIES    - Files to link against to use omniORB
# OMNIORB4_IDL_COMPILER - The idl compiler command (when found)
# OMNIORB4_VERSION      - A string of the form X.Y.Z representing the version
# OMNIORB4_OMNINAMES_COMMAND - the omniNames ORB server command (when found)
# OMNIORB4_DIR          - (optional) Suggested installation directory to search
# 
#
# OMNIORB4_DIR can be used to make it simpler to find the various include
# directories and compiled libraries when omniORB was not installed in the
# usual/well-known directories (e.g. because you made an in tree-source
# compilation or because you installed it in an "unusual" directory).
# Just set OMNIORB4_DIR to point to your specific installation directory.
#
# $ENV{OMNIORB4_DIR} i.e. the environment variable (as opposed to OMNIORB4_DIR
# which is the cmake variable with the same name) is also used to search for
# an omniORB installation.
#
# WARNING: The order of precedence is the following
#    1/ when set OMNIORB4_DIR (the cmake variable)
#    2/ when set $ENV{OMNIORB4_DIR} (the environment variable)
#    3/ the default system pathes !
# This precedence order goes against the usual pratice and default behavior
# of cmake's FIND_* macros. For more on this debate see e.g.
#     http://www.mail-archive.com/kde-buildsystem@kde.org/msg00589.html
#        
#
#### The list of required libraries for omniORB is not so clear. See
# http://www.omniorb-support.com/pipermail/omniorb-list/2005-May/026666.html
# The answer should basically be in
#    - omniORB-4.0.X/README.*
#    - omniORB-4.0.X/readmes/*
# but doesn't seem that obvious to newbies.
# The logic we took is here to impose the detection of omniORB4, omnithread,
# and omniDynamic4 libraries. Optionaly we search for libraries COS4 and
# COSDynamic4 and when present we add them to list of required libraries.
# [the logic is here: the installer of omniORB knows better...]
#

FIND_PATH( OMNIORB4_INCLUDE_DIR
  omniORB4/CORBA.h
  PATHS ${OMNIORB4_DIR}/include $ENV{OMNIORB4_DIR}/include
  NO_DEFAULT_PATH
)
FIND_PATH( OMNIORB4_INCLUDE_DIR omniORB4/CORBA.h )

FIND_LIBRARY( OMNIORB4_LIBRARY_omniORB4
  NAMES omniORB4
  PATHS ${OMNIORB4_DIR}/lib $ENV{OMNIORB4_DIR}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY( OMNIORB4_LIBRARY_omniORB4 NAMES omniORB4 )

FIND_LIBRARY( OMNIORB4_LIBRARY_omnithread
  NAMES omnithread
  PATHS ${OMNIORB4_DIR}/lib $ENV{OMNIORB4_DIR}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY( OMNIORB4_LIBRARY_omnithread NAMES omnithread )

FIND_LIBRARY( OMNIORB4_LIBRARY_omniDynamic4
  NAMES omniDynamic4
  PATHS ${OMNIORB4_DIR}/lib $ENV{OMNIORB4_DIR}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY( OMNIORB4_LIBRARY_omniDynamic4 NAMES omniDynamic4 )

### Optional library
FIND_LIBRARY( OMNIORB4_LIBRARY_COS4
  NAMES COS4
  PATHS ${OMNIORB4_DIR}/lib $ENV{OMNIORB4_DIR}/lib
  NO_DEFAULT_PATH 
)
FIND_LIBRARY( OMNIORB4_LIBRARY_COS4 NAMES COS4 )

### Additional optional library
FIND_LIBRARY( OMNIORB4_LIBRARY_COSDynamic4
  NAMES COSDynamic4
  PATHS ${OMNIORB4_DIR}/lib $ENV{OMNIORB4_DIR}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY( OMNIORB4_LIBRARY_COSDynamic4 NAMES COSDynamic4 )

FIND_PROGRAM( OMNIORB4_IDL_COMPILER
  NAMES omniidl
  PATHS ${OMNIORB4_DIR}/bin $ENV{OMNIORB4_DIR}/bin
  DOC "What is the path where omniidl (the idl compiler) can be found"
  NO_DEFAULT_PATH
)
FIND_PROGRAM( OMNIORB4_IDL_COMPILER
  NAMES omniidl
  DOC "What is the path where omniidl (the idl compiler) can be found"
)

FIND_PROGRAM( OMNIORB4_OMNINAMES_COMMAND
  NAMES omniNames
  PATHS ${OMNIORB4_DIR}/bin $ENV{OMNIORB4_DIR}/bin
  DOC "What is the path where omniNames (the ORB server) can be found"
  NO_DEFAULT_PATH
)
FIND_PROGRAM( OMNIORB4_OMNINAMES_COMMAND
  NAMES omniNames
  DOC "What is the path where omniNames (the ORB server) can be found"
)
 
SET( OMNIORB4_FOUND "NO" )

IF( OMNIORB4_INCLUDE_DIR )
IF( OMNIORB4_LIBRARY_omniORB4 )
IF( OMNIORB4_LIBRARY_omnithread )
IF( OMNIORB4_LIBRARY_omniDynamic4 )
IF( OMNIORB4_IDL_COMPILER )

  SET( OMNIORB4_FOUND "YES" )
  MARK_AS_ADVANCED( OMNIORB4_DIR )
  MARK_AS_ADVANCED( OMNIORB4_INCLUDE_DIR )
  MARK_AS_ADVANCED( OMNIORB4_LIBRARY_omniORB4 )
  MARK_AS_ADVANCED( OMNIORB4_LIBRARY_omnithread )
  MARK_AS_ADVANCED( OMNIORB4_LIBRARY_omniDynamic4 )
  MARK_AS_ADVANCED( OMNIORB4_IDL_COMPILER )
  MARK_AS_ADVANCED( OMNIORB4_OMNINAMES_COMMAND )
  # Note: when linking in static with g++ it looks like omniDynamic4 needs
  # to appear prior to omniORB4 in order to avoid some undefined references
  # to e.g. CORBA::Fixed::Fixed(int).
  SET( OMNIORB4_LIBRARIES
    ${OMNIORB4_LIBRARY_omniDynamic4}
    ${OMNIORB4_LIBRARY_omniORB4}
    ${OMNIORB4_LIBRARY_omnithread}
  )

  # We do not require the existence of the libraries COS4 and COSDynamic4
  # since they exist in some distribution packages (e.g. Fedora Core) and
  # not in others (e.g. Debian). When present we add them to the list of
  # libraries we link against. When absent we hide away (in the advanced
  # parameters) the NOT-FOUND associated variable in order not to deconcert
  # the newbie:
  IF( OMNIORB4_LIBRARY_COS4 )
    SET( OMNIORB4_LIBRARIES
      ${OMNIORB4_LIBRARIES}
      ${OMNIORB4_LIBRARY_COS4}
    )
  ENDIF( OMNIORB4_LIBRARY_COS4 )
  MARK_AS_ADVANCED( OMNIORB4_LIBRARY_COS4 )

  IF( OMNIORB4_LIBRARY_COSDynamic4 )
    SET( OMNIORB4_LIBRARIES
      ${OMNIORB4_LIBRARIES}
      ${OMNIORB4_LIBRARY_COSDynamic4}
    )
  ENDIF( OMNIORB4_LIBRARY_COSDynamic4 )
  MARK_AS_ADVANCED( OMNIORB4_LIBRARY_COSDynamic4 )

  # Optionaly, extract the the version number from the acconfig.h file:
  IF( EXISTS ${OMNIORB4_INCLUDE_DIR}/omniORB4/acconfig.h )
    FILE( READ ${OMNIORB4_INCLUDE_DIR}/omniORB4/acconfig.h OMNIORB_ACCONFIG_H )
    STRING( REGEX MATCH "#define[\t ]+PACKAGE_VERSION[\t ]+\"([0-9]+.[0-9]+.[0-9]+)\"" OMNIORB_ACCONFIG_H "${OMNIORB_ACCONFIG_H}" )
    STRING( REGEX REPLACE ".*\"([0-9]+.[0-9]+.[0-9]+)\".*" "\\1" OMNIORB4_VERSION "${OMNIORB_ACCONFIG_H}" )
  ELSE( EXISTS ${OMNIORB4_INCLUDE_DIR}/omniORB4/acconfig.h )
    SET( OMNIORB4_VERSION "NOT-FOUND" )
  ENDIF( EXISTS ${OMNIORB4_INCLUDE_DIR}/omniORB4/acconfig.h )
  SET( OMNIORB4_VERSION ${OMNIORB4_VERSION}
      CACHE STRING "OmniORB version number." )
  MARK_AS_ADVANCED( OMNIORB4_VERSION )

ENDIF( OMNIORB4_IDL_COMPILER )
ENDIF( OMNIORB4_LIBRARY_omniDynamic4 )
ENDIF( OMNIORB4_LIBRARY_omnithread )
ENDIF( OMNIORB4_LIBRARY_omniORB4 )
ENDIF( OMNIORB4_INCLUDE_DIR )
