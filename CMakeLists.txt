# CMake --------------------------------------------------------------
cmake_minimum_required (VERSION 2.8)
enable_language(C)
enable_language(CXX)
enable_testing()
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(CheckIncludeFiles)
check_include_files(assert.h HAVE_ASSERT_H)
include(CheckIncludeFileCXX)
check_include_file_cxx(cxxabi.h HAVE_CXXABI_H)
add_definitions(-DHAVE_CONFIG_H)

# Project ------------------------------------------------------------
project (codeare CXX C)
set (codeare_VERSION_MAJOR 2)
set (codeare_VERSION_MINOR 0)
set (PROJECT_URL "http://codeare.org/")
set (PROJECT_VERSION
  "${codeare_VERSION_MAJOR}.${codeare_VERSION_MINOR}")
set (PACKAGE_NAME "codeare")
set (PACKAGE_VERSION "${PROJECT_VERSION}")
set (top_srcdir "${PROJECT_SOURCE_DIR}")

# Git ----------------------------------------------------------------
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# OS type ------------------------------------------------------------
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set (Boost_USE_STATIC_LIBS ON)
  set(WINDOWS TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(LINUX TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX TRUE)
endif()

# Architecture -------------------------------------------------------
include (VcMacros)
include (OptimizeForArchitecture)

vc_determine_compiler()

if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "(i686|x86|AMD64|amd64)")
  set(Vc_X86 TRUE)
elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "(arm|aarch32|aarch64)")
  message(WARNING "No optimized implementation of the Vc types available for ${CMAKE_SYSTEM_PROCESSOR}")
  set(Vc_ARM TRUE)
else()
  message(WARNING "No optimized implementation of the Vc types available for ${CMAKE_SYSTEM_PROCESSOR}")
endif()

option(USE_CCACHE "If enabled, ccache will be used (if it exists on the system) to speed up recompiles." OFF)
if(USE_CCACHE)
  find_program(CCACHE_COMMAND ccache)
  if(CCACHE_COMMAND)
    mark_as_advanced(CCACHE_COMMAND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_COMMAND}")
  endif()
endif()

# TODO: check that 'decltype' compiles
# TODO: check that 'constexpr' compiles
if(NOT Vc_COMPILER_IS_MSVC) # MSVC doesn't provide a switch to turn C++11 on/off AFAIK
  AddCompilerFlag("-std=c++14" CXX_RESULT _ok CXX_FLAGS CMAKE_CXX_FLAGS)
  if(NOT _ok)
    AddCompilerFlag("-std=c++1y" CXX_RESULT _ok CXX_FLAGS CMAKE_CXX_FLAGS)
    if(NOT _ok)
      AddCompilerFlag("-std=c++11" CXX_RESULT _ok CXX_FLAGS CMAKE_CXX_FLAGS)
      if(NOT _ok)
        AddCompilerFlag("-std=c++0x" CXX_RESULT _ok CXX_FLAGS CMAKE_CXX_FLAGS)
        if(NOT _ok)
          message(FATAL_ERROR "Vc 1.x requires C++11, better even C++14. It seems this is not available. If this was incorrectly determined please notify vc-devel@compeng.uni-frankfurt.de")
        endif()
      endif()
    endif()
  endif()
elseif(Vc_MSVC_VERSION LESS 180021114)
  message(FATAL_ERROR "Vc 1.x requires C++11 support. This requires at least Visual Studio 2013 with the Nov 2013 CTP.")
endif()

if(Vc_COMPILER_IS_GCC)
  if(Vc_GCC_VERSION VERSION_GREATER "5.0.0" AND Vc_GCC_VERSION VERSION_LESS "6.0.0")
    UserWarning("GCC 5 goes into an endless loop comiling example_scaling_scalar. Therefore, this target is disabled.")
    list(APPEND disabled_targets
      example_scaling_scalar
      )
  endif()
elseif(Vc_COMPILER_IS_MSVC)
  if(MSVC_VERSION LESS 1700)
    # MSVC before 2012 has a broken std::vector::resize implementation. STL + Vc code will probably not compile.
    # UserWarning in VcMacros.cmake
    list(APPEND disabled_targets
      stlcontainer_sse
      stlcontainer_avx
      )
  endif()
  # Disable warning "C++ exception specification ignored except to indicate a function is not __declspec(nothrow)"
  # MSVC emits the warning for the _UnitTest_Compare desctructor which needs the throw declaration so that it doesn't std::terminate
  AddCompilerFlag("/wd4290")
endif()

vc_set_preferred_compiler_flags(BUILDTYPE_FLAGS)

add_definitions(${Vc_DEFINITIONS})
add_compile_options(${Vc_COMPILE_FLAGS})

if(Vc_COMPILER_IS_INTEL)
  # per default icc is not IEEE compliant, but we need that for verification
  AddCompilerFlag("-fp-model source")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
  message(STATUS "WARNING! It seems you are compiling without optimization. Please set CMAKE_BUILD_TYPE.")
endif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")


# C++ flags ----------------------------------------------------------
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} -Wno-psabi -DTIXML_USE_STL -fPIC -DHAVE_CXXABI_H")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  if (${MACOSX})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL -DHAVE_CXXABI_H")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL -fPIC -DHAVE_CXXABI_H")
  endif()
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_STL -DHAVE_CXXABI_H")
elseif (${MSVC})
  set(CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} /DTIXML_USE_STL /EHsc /Ox /nologo /wd4267 /wd4244 /wd4190 /wd4996 /wd4251 /wd4305 /LD /MT /DEXP_STL") 
  set (CMAKE_LINKER_FLAGS 
    "${CMAKE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_SHARED_LINKER_FLAGS 
    "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set (CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  if (${AVX_FOUND})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  elseif (${SSE2_FOUND})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
  endif()
endif()

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("emmintrin.h" HAVE_EMMINTRIN_H)

# C++ 0x/11 ----------------------------------------------------------
include (CheckCXX11Features)
if (${HAS_CXX11_SHARED_PTR}) 
  add_definitions (-DHAVE_CXX11_SHARED_PTR)
endif ()
if (${HAS_CXX11_RVALUE_REFERENCES}) 
  add_definitions (-DHAVE_CXX11_RVALUE_REFERENCES)
endif ()
if (${HAS_CXX11_MUTEX}) 
  add_definitions (-DHAVE_CXX11_MUTEX)
endif ()
if (${HAS_CXX11_NOEXCEPT}) 
  add_definitions (-DHAVE_CXX11_NOEXCEPT)
endif ()
if (${HAS_CXX11_TUPLE}) 
  add_definitions (-DHAVE_CXX11_TUPLE)
endif ()
if (${HAS_CXX11_CONDITIONAL}) 
  add_definitions (-DHAVE_CXX11_CONDITIONAL)
endif ()
if (${HAS_CXX11_THREAD})
  add_definitions (-DHAVE_CXX11_THREAD)
endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")

# OpenMP  -------------------------------------------------------
find_package (OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Core ---------------------------------------------------------------
find_package(Boost COMPONENTS chrono timer system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBoost_MAJOR_VERSION=${Boost_MAJOR_VERSION})
add_definitions(-DBoost_MINOR_VERSION=${Boost_MINOR_VERSION})
link_directories (${BOOST_LIBRARYDIR})
#find_package(OpenSSL REQUIRED)
#include_directories(${OPENSSL_INCLUDE_DIR})
#set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OPENSSL_LIBRARIES}")


# MATLAB --------------------------------------------------------------
if (BUILD_MATLAB_SUPPORT) 
  find_package(Matlab COMPONENTS MX_LIBRARY MAT_LIBRARY)
  if(Matlab_ROOT_DIR)
    include_directories (${Matlab_INCLUDE_DIRS})
    add_definitions (-DHAVE_MAT_H)
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
  include_directories (${OMNIORB4_INCLUDE_DIR})
else ()
  message ("OmniORB4 library not found. Networking infrastructure will not bebuilt.")
endif ()

# Fourier transform --------------------------------------------------
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
find_package(NFFT3 COMPONENTS double threads)
if (NFFT3_FOUND)
  add_definitions(-DHAVE_NFFT3)
  set (NFFT3_LIBRARIES_TMP ${NFFT3_LIBRARIES})
  set (NFFT3_INCLUDE_DIR_TMP ${NFFT3_INCLUDE_DIR})
  unset (NFFT3_LIBRARIES)
  unset (NFFT3_INCLUDE_DIR)
  unset (NFFT3_FOUND)
  find_package(NFFT3 COMPONENTS single double threads)
  if (NFFT3_FOUND MATCHES FALSE)
    set (NFFT3_FOUND TRUE)
    set (NFFT3_LIBRARIES ${NFFT3_LIBRARIES_TMP})
    set (NFFT3_INCLUDE_DIR ${NFFT3_INCLUDE_DIR_TMP})
  else ()
    add_definitions(-DHAVE_NFFT3F)
  endif()
  include_directories (${NFFT3_INCLUDE_DIR})
  include (CheckTypeSize)
  include(CheckStructHasMember)
  check_struct_has_member(solver_plan_complex flags nfft3.h SOLVER_PLAN_COMPLEX_HAS_MEMBER_FLAGS)
  if (${SOLVER_PLAN_COMPLEX_HAS_MEMBER_FLAGS})
    add_definitions(-DNFFT_HAS_MEMBER_FLAGS)
  endif()
  check_type_size(nfft_mv_plan_complex NFFT_MV_PLAN_COMPLEX_SIZE)
  if (${NFFT_MV_PLAN_COMPLEX_SIZE})
    add_definitions(-DUSE_NFFT_32_NAMING)
  endif()
endif()
if (NOT NFFT3_INCLUDE_DIR)
   set(NFFT3_INCLUDE_DIR "")
endif()
# Interpolations -----------------------------------------------------
find_package (PNG)
find_package (TIFF)
if (${PNG_FOUND} AND ${TIFF_FOUND})
  find_package(ITK 4.4 COMPONENTS ITKCommon ITKIOImageBase ITKIOGDCM)
  if (${ITK_FOUND})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DITK_VERSION_MAJOR=${ITK_VERSION_MAJOR}")
    include_directories (${ITK_INCLUDE_DIRS})
    message (STATUS "Found ITK: ${ITK_USE_FILE} ${ITK_LIBRARY_DIRS}")
    add_definitions (-DVCL_CAN_STATIC_CONST_INIT_FLOAT=0)
  else()
    message (WARNING "Couldn't find ITK. Image resampling and DICOM IO will not be built.") 
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
find_package(HDF5 "1.8.5" REQUIRED COMPONENTS C CXX)
include_directories (${HDF5_INCLUDE_DIRS})
link_directories (${HDF5_LIBRARY_DIRS})
string(REPLACE "." ";" HDF5_VERSION_LIST ${HDF5_VERSION})
list(GET HDF5_VERSION_LIST 0 HDF5_VERSION_MAJ)
list(GET HDF5_VERSION_LIST 1 HDF5_VERSION_MIN)
list(GET HDF5_VERSION_LIST 2 HDF5_VERSION_REL)
if ((HDF5_VERSION_MAJ GREATER 0) AND (HDF5_VERSION_MIN GREATER 7) AND (HDF5_VERSION_REL GREATER 11))
  add_definitions(-DHDF5_HAS_LOCATION)
endif() 
if (${MSVC})
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} hdf5.lib hdf5_cpp.lib")
  set(ZLIB_LIBRARY "${HDF5_DIR}/../../lib/zlib.lib")
  set(ZLIB_INCLUDE_DIR "${HDF5_DIR}/../../include")
  set(ZLIB_FOUND "1")
else()
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${HDF5_LIBRARIES}")
  find_package(ZLIB REQUIRED)
  list (APPEND HDF5_LIBRARIES "z")
endif()

find_package(Ismrmrd)
if (${ISMRMRD_FOUND})
  include_directories(${ISMRMRD_INCLUDE_DIR})
  link_directories(${ISMRMRD_INCLUDE_DIR})
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

add_subdirectory(doc)
add_subdirectory(src)
