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

add_executable (t_codeare t_codeare.cpp)
add_test (cod t_codeare)

add_executable (t_hdf5 t_hdf5.cpp)
add_test (hdf5 t_hdf5)

if (${MATLAB_FOUND})
	add_executable (t_matlab t_matlab.cpp)
	target_link_libraries (t_matlab ${MATLAB_LIBRARIES})
	add_test (matlab t_matlab)
endif()

if (${NIFTI_FOUND})
	add_executable (t_nifti t_nifti.cpp)
	target_link_libraries (t_nifti ${NIFTI_NIFTIIO_LIBRARY} ${NIFTI_ZNZ_LIBRARY} ${ZLIB_LIBRARIES})
	add_test (nifti t_nifti)
endif() 

