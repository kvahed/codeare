cmake_minimum_required (VERSION 2.6)

project (codeare)
set (codeare_VERSION_MAJOR 1)
set (codeare_VERSION_MINOR 2)

set(PROJECT_URL "http://codeare.org/")
set(PROJECT_VERSION
"${codeare_VERSION_MAJOR}.${codeare_VERSION_MINOR}")

enable_testing()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

include(ConfigureChecks.cmake)
add_definitions(-DHAVE_CONFIG_H)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox /nologo /MT")
endif()

# Architecture
include (VcMacros)
include (OptimizeForArchitecture)
find_package (OpenMP)

# Core
find_package(Boost COMPONENTS thread system REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(ZLIB REQUIRED)

# Linear algebra
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Fourier transform
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
find_package(NFFT)

# Interpolations
find_package (PNG)
find_package (TIFF)
if (${PNG_FOUND} AND ${TIFF_FOUND})
  find_package(ITK)
endif()
find_package(GSL)

# File IO 
find_package(HDF5 COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_LINK_FLAGS "-L${HDF5_LIBRARY_DIRS} -lhdf5 -lhdf5_cpp")
find_package(Ismrmrd REQUIRED)
find_package(Matlab)
if (${ZLIB_FOUND}) 
  find_package(NIFTI) 
endif()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
