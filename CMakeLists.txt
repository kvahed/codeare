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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox
  /nologo /MT /wd4267 /wd4244 /wd4190 /wd4996 /LD /MD /DEXP_STL")
endif()

# Architecture ------------------------------------------------------------------
include (VcMacros)
include (OptimizeForArchitecture)
find_package (OpenMP)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}") 
endif()

# Core --------------------------------------------------------------------------
find_package(Boost COMPONENTS thread system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OPENSSL_LIBRARIES}")


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
find_package (PNG)
find_package (TIFF)
if (${PNG_FOUND} AND ${TIFF_FOUND})
  find_package(ITK)
endif()
find_package(GSL)

# File IO -----------------------------------------------------------------------
find_package(HDF5 COMPONENTS CXX REQUIRED)
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

find_package(Matlab)
if(${MATLAB_FOUND})
  include_directories (${MATLAB_INCLUDE_DIRS})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    link_directories ("${MATLAB_ROOT_DIR}/extern/lib/win64/microsoft")
  endif()

endif()

if (${ZLIB_FOUND}) 
  find_package(NIFTI) 
  if (${NIFTI_INCLUDE_DIR} MATCHES "")
    set(NIFTI_INCLUDE_DIR "")
  endif()
endif()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
