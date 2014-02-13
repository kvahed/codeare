AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/mri -I$(top_srcdir)/src/matrix/arithmetic -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix/io -I$(top_srcdir)/src/core

TESTS = t_phase_combine
check_PROGRAMS = $(TESTS)

t_phase_combine_SOURCES = t_phase_combine.cpp
t_phase_combine_LDFLAGS =  $(BLAS) $(LAPACK)

