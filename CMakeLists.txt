cmake_minimum_required (VERSION 2.6)

project (codeare)
set (codeare_VERSION_MAJOR 1)
set (codeare_VERSION_MINOR 2)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(warnings "")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(warnings "")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox /nologo /MT")
endif()

include (VcMacros)
include (OptimizeForArchitecture)

find_package(Boost COMPONENTS thread system REQUIRED)
find_package(HDF5 1.8 COMPONENTS C CXX REQUIRED)
find_package(BLAS REQUIRED)
find_package(FFTW3 COMPONENTS double single threads)
find_package(OpenSSL) 

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
