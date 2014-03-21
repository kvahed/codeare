# Project ------------------------------------------------------------
project (codeare CXX C)
set (codeare_VERSION_MAJOR 1)
set (codeare_VERSION_MINOR 2)
set (PROJECT_URL "http://codeare.org/")
set (PROJECT_VERSION
  "${codeare_VERSION_MAJOR}.${codeare_VERSION_MINOR}")
set (PACKAGE_NAME "codeare")
set (PACKAGE_VERSION "${PROJECT_VERSION}")
set (top_srcdir "${PROJECT_SOURCE_DIR}")

# CMake --------------------------------------------------------------
cmake_minimum_required (VERSION 2.8)
enable_testing()
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(CheckIncludeFiles)
check_include_files(assert.h HAVE_ASSERT_H)
add_definitions(-DHAVE_CONFIG_H)

# Git ----------------------------------------------------------------
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# OS type ------------------------------------------------------------
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(WINDOWS TRUE)
  set (Boost_USE_STATIC_LIBS ON)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(LINUX TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX TRUE)
  set (Boost_USE_STATIC_LIBS ON)
endif()

# C++ flags ----------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} -Wno-psabi -DTIXML_USE_STL -fPIC -Ofast")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox /nologo /wd4267 /wd4244 /wd4190 /wd4996 /LD /MD /DEXP_STL") 
  set (CMAKE_LINKER_FLAGS 
    "${CMAKE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_SHARED_LINKER_FLAGS 
    "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
endif()

# C++ 0x/11 ----------------------------------------------------------
include (CheckCXX11Features)
if (${HAS_CXX11_SHARED_PTR}) 
  add_definitions (-DHAVE_CXX11_SHARED_PTR)
endif ()
#if (${HAS_CXX11_THREAD})  IMMATURE IMPL?
#  add_definitions (-DHAVE_CXX11_THREAD)
#endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")

# Architecture -------------------------------------------------------
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

# Core ---------------------------------------------------------------
find_package(Boost COMPONENTS thread chrono timer system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories (${BOOST_LIBRARYDIR})
find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OPENSSL_LIBRARIES}")


# MATLAB --------------------------------------------------------------
find_package(Matlab)
if(${MATLAB_FOUND})
  include_directories (${MATLAB_INCLUDE_DIRS})
  add_definitions (-DHAVE_MAT_H)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    link_directories ("${MATLAB_ROOT_DIR}/extern/lib/win64/microsoft")
  else()
    link_directories ("${MATLAB_LIBRARY_DIRS}")
  endif()
endif()

# Linear algebra -----------------------------------------------------
find_package(LAPACK REQUIRED)
if (${WINDOWS})
   link_directories (${INTEL_COMPILER_DIR}/lib/intel64)
endif ()

# OmniORB ------------------------------------------------------------
find_package (OmniORB)
if (${OMNIORB4_FOUND})
  message ("OmniORB4 library not found. Networking infrastructure will not bebuilt.")
  include_directories (${OMNIORB4_INCLUDE_DIR})
endif ()

# Fourier transform --------------------------------------------------
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
find_package(NFFT)
if(${NFFT_INCLUDE_DIR} MATCHES "")
  set (NFFT_INCLUDE_DIR "")
endif()

# Interpolations -----------------------------------------------------
find_package (PNG)
find_package (TIFF)
if (${PNG_FOUND} AND ${TIFF_FOUND})
  find_package(ITK)
  if (${ITK_FOUND})
    include_directories (${ITK_INCLUDE_DIRS})
    link_directories (${ITK_LIBRARY_DIRS})
  endif()
endif()

find_package(GSL)
if (${GSL_FOUND})
  include_directories (${GSL_INCLUDE_DIR})
  add_definitions (-DHAVE_GSL)
endif()

# OpenCL -------------------------------------------------------------
find_package (OpenCL)
if (OPENCL_FOUND)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

# MPI ----------------------------------------------------------------
find_package (MPI)
if (${MPI_FOUND})
  message (STATUS "Found MPI.")
  include_directories(${MPI_INCLUDE_DIRS})
  #find_package (Scalapack)
  #if (${SCALAPACK_FOUND})
  #  message (STATUS "Found Scalapack.")
  #endif()
endif()

# File IO ------------------------------------------------------------
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
  message (WARNING 
    "Couldn't find ISMRMRD installation. Support will not be built in.") 
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
add_subdirectory(doc)
