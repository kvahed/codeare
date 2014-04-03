########################################
# BEGIN_COPYRIGHT
#
# This file is part of SciDB.
# Copyright (C) 2008-2013 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the AFFERO GNU General Public License as published by
# the Free Software Foundation.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the AFFERO GNU General Public License for the complete license terms.
#
# You should have received a copy of the AFFERO GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
#
# END_COPYRIGHT
########################################

# SCALAPACK_FOUND
# SCALAPACK_LIBRARIES

set(CMAKE_SYSTEM_PREFIX_PATH_BACKUP ${CMAKE_SYSTEM_PREFIX_PATH})

#
# setup for the find_library
#
if(${DISTRO_NAME_VER} MATCHES "RedHat-[0-9][.][0-9]")
  # RedHat 5.4
  set(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH} "/usr/lib64/openmpi/1.4-gcc")
  set(SCALAPACK_LIB_NAME scalapack) # name for search
endif()

if(${DISTRO_NAME_VER} MATCHES "Fedora")
  # Fedora 11..17
  set(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH} "/usr/lib64/openmpi")
  set(SCALAPACK_LIB_NAME scalapack) # name for search
  set(BLACS_REQUIRED 1)
  set(BLACS_LIB_NAME mpiblacs) # name for search
  set(BLACS_C_INIT_REQUIRED 1)
  set(BLACS_C_INIT_LIB_NAME mpiblacsCinit) # name for search
endif()

if(${DISTRO_NAME_VER} MATCHES "Fedora-11")
  # Fedora 11
  set(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH} "/usr/lib64/scalapack-openmpi")
endif()

if(${DISTRO_NAME_VER} MATCHES "Ubuntu")
  set(SCALAPACK_LIB_NAME scalapack-openmpi) # name for search
  set(BLACS_REQUIRED 1)
  set(BLACS_LIB_NAME blacs-openmpi) # name for search
endif()

#
# do the find_library
#
if(NOT DISABLE_SCALAPACK)

  find_library(SCALAPACK_LIB NAMES ${SCALAPACK_LIB_NAME})

  if(BLACS_REQUIRED)
    find_library(BLACS_LIB NAMES ${BLACS_LIB_NAME})
  endif()

  if(BLACS_C_INIT_REQUIRED)
    find_library(BLACS_C_INIT_LIB NAMES ${BLACS_C_INIT_LIB_NAME})
  endif()


  if(${DISTRO_NAME_VER} MATCHES "Ubuntu")
    # Ubuntu
    if ("${BLACS_LIB}" STREQUAL "BLACS_LIB-NOTFOUND")
      # Ubuntu 11.04 has only libblacs-openmpi.so.1 which find_library does not match
      find_file(BLACS_LIB libblacs-openmpi.so.1 PATHS /usr/lib)
    endif ("${BLACS_LIB}" STREQUAL "BLACS_LIB-NOTFOUND")

    #
    # fix repairable cases of failed find_library
    #
    if ("${SCALAPACK_LIB}" STREQUAL "SCALAPACK_LIB-NOTFOUND")
      # Ubuntu 11.04 has only /usr/lib/libscalapack-openmpi.so.1 only, find_library() misses it
      find_file(SCALAPACK_LIB libscalapack-openmpi.so.1 PATHS /usr/lib)
    endif ("${SCALAPACK_LIB}" STREQUAL "SCALAPACK_LIB-NOTFOUND")
  endif()

  set(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH_BACKUP})

  #
  # check for completeness and set outputs
  #
  if (NOT ("${SCALAPACK_LIB}" STREQUAL "SCALAPACK_LIB-NOTFOUND"))
     set (SCALAPACK_LIB_FOUND 1)
  endif(NOT ("${SCALAPACK_LIB}" STREQUAL "SCALAPACK_LIB-NOTFOUND"))

  if(BLACS_REQUIRED)
    if(NOT "${BLACS_LIB}" STREQUAL "BLACS_LIB-NOTFOUND")
      set(BLACS_LIB_FOUND 1)
    endif (NOT "${BLACS_LIB}" STREQUAL "BLACS_LIB-NOTFOUND")
  endif(BLACS_REQUIRED)

  if(BLACS_C_INIT_REQUIRED)
    if(NOT "${BLACS_C_INIT_LIB}" STREQUAL "BLACS_C_INIT_LIB-NOTFOUND")
      set(BLACS_C_INIT_LIB_FOUND 1)
    endif(NOT "${BLACS_C_INIT_LIB}" STREQUAL "BLACS_C_INIT_LIB-NOTFOUND")
  endif(BLACS_C_INIT_REQUIRED)

  #MESSAGE(STATUS "SCALAPACK_LIB_FOUND=${SCALAPACK_LIB_FOUND}")
  #MESSAGE(STATUS "SCALAPACK_LIBD=${SCALAPACK_LIB}")

  #MESSAGE(STATUS "BLACS_REQUIRED=${BLACS_REQUIRED}")
  #MESSAGE(STATUS "BLACS_LIB_FOUND=${BLACS_LIB_FOUND}")
  #MESSAGE(STATUS "BLACS_LIB=${BLACS_LIB}")

  #MESSAGE(STATUS "BLACS_C_INIT_REQUIRED=${BLACS_REQUIRED}")
  #MESSAGE(STATUS "BLACS_C_INIT_LIB_FOUND=${BLACS_C_INIT_LIB_FOUND}")
  #MESSAGE(STATUS "BLACS_C_INIT_LIB=${BLACS_C_INIT_LIB}")

  if(SCALAPACK_LIB_FOUND AND
     ((NOT BLACS_REQUIRED) OR BLACS_LIB_FOUND) AND
     ((NOT BLACS_C_INIT_REQUIRED) OR BLACS_C_INIT_LIB_FOUND))

    # ScaLAPACK
    set(SCALAPACK_LIBRARIES "${SCALAPACK_LIB}")

    # BLACS
    if(BLACS_REQUIRED)
      set(SCALAPACK_LIBRARIES "${BLACS_LIB};${SCALAPACK_LIBRARIES}")
    endif(BLACS_REQUIRED)

    # BLACS cinit
    if(BLACS_C_INIT_REQUIRED)
      set(SCALAPACK_LIBRARIES "${BLACS_C_INIT_LIB};${SCALAPACK_LIBRARIES}")
    endif(BLACS_C_INIT_REQUIRED)

    set(SCALAPACK_FOUND 1)

  endif(SCALAPACK_LIB_FOUND AND
    ((NOT BLACS_REQUIRED) OR BLACS_LIB_FOUND) AND
    ((NOT BLACS_C_INIT_REQUIRED) OR BLACS_C_INIT_LIB_FOUND))

endif(NOT DISABLE_SCALAPACK)
