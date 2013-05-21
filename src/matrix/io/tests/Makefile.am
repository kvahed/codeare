AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/io 
# -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix/ft -I$(top_srcdir)/src/misc

TESTS = t_hdf5 t_matlab
check_PROGRAMS = $(TESTS)

t_hdf5_SOURCES = $(COMMON) t_hdf5.cpp
t_matlab_SOURCES = $(COMMON) t_matlab.cpp


