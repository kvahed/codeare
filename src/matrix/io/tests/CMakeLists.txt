include_directories(
  ${PROJECT_SOURCE_DIR}/src/core
  ${PROJECT_SOURCE_DIR}/src/matrix
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/matrix/simd
  ${HDF5_INCLUDE_DIRS}
  ${MATLAB_INCLUDE_DIRS}
  if (${NIFTI_FOUND})
    ${NIFTI_INCLUDE_DIR}
  endif()
  }
)

include (TestMacro)

if (${LINUX})
  set (LD_ENV "LD_LIBRARY_PATH=$LD_LIBRARY_PATH")
elseif (${MACOSX})
  set (LD_ENV "DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH")
elseif (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
  list (APPEND LD_ENV ${HDF5_TOOLS_DIR})
  list (APPEND LD_ENV Release)
endif()

add_executable (t_codeare t_codeare.cpp)
add_test (cod t_codeare)

add_executable (t_vxfile t_vxfile.cpp)
add_test (vx t_vxfile)
target_link_libraries (t_vxfile ${OPENSSL_LIBRARIES} boost_timer boost_system core)

# if (${ITK_FOUND})
#   add_executable (t_dicom t_dicom.cpp)
#   add_test (dicom t_dicom)
#   target_link_libraries (t_dicom ${OPENSSL_LIBRARIES} itksys ITKCommon ITKIOImageBase ITKIOGDCM)
# endif ()

add_executable (t_hdf5 t_hdf5.cpp)
if (!${WINDOWS})
  target_link_libraries (t_hdf5 ${HDF5_LIBRARIES} ${OPENSSL_LIBRARIES} core)
else()
  target_link_libraries (t_hdf5 ${HDF5_LIBRARIES} ${OPENSSL_LIBRARIES} core)
endif()
set (TEST_CALL t_hdf5)  
MP_TESTS ("hdf5" "${TEST_CALL}")
set_tests_properties (hdf5 PROPERTIES ENVIRONMENT
  "HDF5_DISABLE_VERSION_CHECK=2") 

if(Matlab_ROOT_DIR)
  add_executable (t_matlab t_matlab.cpp)
  target_link_libraries (t_matlab ${Matlab_MX_LIBRARY} ${Matlab_MAT_LIBRARY})
  
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${MATLAB_ROOT}/bin/win64") 
  else ()
    set (LD_ENV "${LD_ENV}:${MATLAB_LIBRARY_DIRS}")
  endif ()

  set (TEST_CALL t_matlab -b 2)  
  MP_TESTS ("matlab" "${TEST_CALL}")
  
endif()

if (${NIFTI_FOUND})
	add_executable (t_nifti t_nifti.cpp)
	target_link_libraries (t_nifti ${NIFTI_NIFTIIO_LIBRARY}
	${NIFTI_ZNZ_LIBRARY} ${ZLIB_LIBRARIES}) 
	add_test (nifti t_nifti)
endif() 

