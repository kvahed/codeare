# - Try to find NFFT3.
# Usage: find_package(NFFT3 [COMPONENTS [single double long-double threads]])
#
# Variables used by this module:
#  NFFT3_ROOT_DIR             - NFFT3 root directory
# Variables defined by this module:
#  NFFT3_FOUND                - system has NFFT3
#  NFFT3_INCLUDE_DIR          - the NFFT3 include directory (cached)
#  NFFT3_INCLUDE_DIRS         - the NFFT3 include directories
#                               (identical to NFFT3_INCLUDE_DIR)
#  NFFT3[FL]?_LIBRARY         - the NFFT3 library - double, single(F), 
#                               long-double(L) precision (cached)
#  NFFT3[FL]?_THREADS_LIBRARY - the threaded NFFT3 library - double, single(F), 
#                               long-double(L) precision (cached)
#  NFFT3_LIBRARIES            - list of all NFFT3 libraries found

# Copyright (C) 2009-2010
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id: FindNFFT3.cmake 15918 2010-06-25 11:12:42Z kvahed $

# Use double precision by default.
if(NFFT3_FIND_COMPONENTS MATCHES "^$")
  set(_components double)
else()
  set(_components ${NFFT3_FIND_COMPONENTS})
endif()

# Loop over each component.
set(_libraries)
foreach(_comp ${_components})
  if(_comp STREQUAL "single")
    list(APPEND _libraries nfft3f)
  elseif(_comp STREQUAL "double")
    list(APPEND _libraries nfft3)
  elseif(_comp STREQUAL "long-double")
    list(APPEND _libraries nfft3l)
  elseif(_comp STREQUAL "threads")
    set(_use_threads ON)
  else(_comp STREQUAL "single")
    message(FATAL_ERROR "FindNFFT3: unknown component `${_comp}' specified. "
      "Valid components are `single', `double', `long-double', and `threads'.")
  endif(_comp STREQUAL "single")
endforeach(_comp ${_components})

# If using threads, we need to link against threaded libraries as well.
if(_use_threads)
  set(_thread_libs)
  foreach(_lib ${_libraries})
    list(APPEND _thread_libs ${_lib}_threads)
  endforeach(_lib ${_libraries})
  #  set(_libraries ${_thread_libs} ${_libraries})
  set(_libraries ${_thread_libs})
endif(_use_threads)

# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set(_check_list)

# Search for all requested libraries.
foreach(_lib ${_libraries})
  string(TOUPPER ${_lib} _LIB)
  find_library(${_LIB}_LIBRARY ${_lib}
    HINTS ${NFFT3_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(${_LIB}_LIBRARY)
  list(APPEND NFFT3_LIBRARIES ${${_LIB}_LIBRARY})
  list(APPEND _check_list ${_LIB}_LIBRARY)
endforeach(_lib ${_libraries})

# Search for the header file.
find_path(NFFT3_INCLUDE_DIR nfft3.h 
  HINTS ${NFFT3_ROOT_DIR} PATH_SUFFIXES include)
mark_as_advanced(NFFT3_INCLUDE_DIR)
list(APPEND _check_list NFFT3_INCLUDE_DIR)

# Handle the QUIETLY and REQUIRED arguments and set NFFT_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NFFT3 DEFAULT_MSG ${_check_list})
