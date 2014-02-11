cmake_minimum_required (VERSION 2.6)

project (codeare)
set (codeare_VERSION_MAJOR 1)
set (codeare_VERSION_MINOR 2)

enable_testing()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

include(ConfigureChecks.cmake)
add_definitions(-DHAVE_CONFIG_H)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

if ($CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(warnings "")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox /nologo /MT")
endif()

include (VcMacros)
include (OptimizeForArchitecture)

# Core
find_package(Boost COMPONENTS thread system REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(ZLIB REQUIRED)

# Linear algebra
find_package(BLAS REQUIRED)
find_package(CBLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Fourier transform
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
#find_package(NFFT3 COMPONENTS single OPTIONAL)

# Interpolations
find_package(ITK)
find_package(GSL)

# File IO 
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(Ismrmrd REQUIRED)
find_package(Matlab)
if (${ZLIB_FOUND}) 
	find_package(NIFTI) 
endif()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
