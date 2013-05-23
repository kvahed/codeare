AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/io 

TESTS = t_hdf5 
if USE_MATLAB
	TESTS += t_matlab
endif 

check_PROGRAMS = $(TESTS)

t_hdf5_SOURCES = $(COMMON) t_hdf5.cpp
if USE_MATLAB
	t_matlab_SOURCES = $(COMMON) t_matlab.cpp
endif 


