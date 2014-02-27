# Project -----------------------------------------------------------------------
project (codeare)
set (codeare_VERSION_MAJOR 1)
set (codeare_VERSION_MINOR 2)
set (PROJECT_URL "http://codeare.org/")
set (PROJECT_VERSION
  "${codeare_VERSION_MAJOR}.${codeare_VERSION_MINOR}")

# CMake -------------------------------------------------------------------------
cmake_minimum_required (VERSION 2.8)
enable_testing()
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(ConfigureChecks.cmake)
add_definitions(-DHAVE_CONFIG_H)

# Git ---------------------------------------------------------------------------
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Configuration for C++ compile time --------------------------------------------
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

# OS type -----------------------------------------------------------------------
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(WINDOWS TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(LINUX TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX TRUE)
endif()

# C++ flags ---------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi -DTIXML_USE_STL -fPIC") 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox
  /nologo /MT /wd4267 /wd4244 /wd4190 /wd4996 /LD /MD /DEXP_STL") 
  set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
endif()

# Architecture ------------------------------------------------------------------
include (VcMacros)
include (OptimizeForArchitecture)
find_package (OpenMP)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
# Core --------------------------------------------------------------------------
find_package(Boost COMPONENTS thread system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories (${BOOST_LIBRARYDIR})
find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OPENSSL_LIBRARIES}")

# OmniORB -----------------------------------------------------------------------
find_package (OmniORB)
if (!{$OMNIORB4_FOUND})
  message (WARNING "OmniORB4 library not found. Networking infrastructure will not bebuilt.")
endif ()

# Linear algebra ----------------------------------------------------------------
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Fourier transform -------------------------------------------------------------
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
find_package(NFFT)
if(${NFFT_INCLUDE_DIR} MATCHES "")
  set (NFFT_INCLUDE_DIR "")
endif()

# Interpolations ----------------------------------------------------------------
if (${WINDOWS})
  find_package (PNG)
  find_package (TIFF)
endif()
find_package(ITK)
find_package(GSL)

# File IO -----------------------------------------------------------------------
find_package(HDF5 REQUIRED)
message (STATUS ${HDF5_LIBRARIES})
include_directories (${HDF5_INCLUDE_DIRS})
link_directories (${HDF5_LIBRARY_DIRS})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} hdf5.lib hdf5_cpp.lib")
  set(ZLIB_LIBRARY "${HDF5_DIR}/../../lib/zlib.lib")
  set(ZLIB_INCLUDE_DIR "${HDF5_DIR}/../../include")
  set(ZLIB_FOUND "1")
else()
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -lhdf5 -lhdf5_cpp")
  find_package(ZLIB REQUIRED)
endif()

find_package(Ismrmrd)
if (${ISMRMRD_FOUND})
  include_directories(${ISMRMRD_INCLUDE_DIR})
  link_directories(${ISMRMRD_INCLUDE_DIR})
else ()
  message (WARNING "Couldn't find ISMRMRD installation. Support will not be built in.") 
endif()

find_package(Matlab)
if(${MATLAB_FOUND})
  include_directories (${MATLAB_INCLUDE_DIRS})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    link_directories ("${MATLAB_ROOT_DIR}/extern/lib/win64/microsoft")
  else()
    link_directories ("${MATLAB_LIBRARY_DIRS}")
  endif()
endif()

if (${ZLIB_FOUND}) 
  find_package(NIFTI) 
  if (${NIFTI_FOUND})
    include_directories (${NIFTI_INCLUDE_DIR})
  else()
    set(NIFTI_INCLUDE_DIR "")   
  endif()
endif()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
