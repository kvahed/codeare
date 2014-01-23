AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/matrix/io 

TESTS = t_codeare t_hdf5 t_iocontext
if USE_NIFTI
TESTS += t_nifti
endif
if USE_MATLAB
TESTS += t_matlab
endif 

check_PROGRAMS = $(TESTS)

t_codeare_SOURCES = $(COMMON) t_codeare.cpp
t_hdf5_SOURCES = $(COMMON) t_hdf5.cpp
t_iocontext_SOURCES = $(COMMON) t_iocontext.cpp
if USE_NIFTI
t_nifti_SOURCES = $(COMMON) t_nifti.cpp
endif
if USE_MATLAB
t_matlab_SOURCES = $(COMMON) t_matlab.cpp
endif 

clean-local:
	rm -f test.mat test.h5 test.cod test.nii



