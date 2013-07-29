AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/ft -I$(top_srcdir)/src/matrix/arithmetic -I$(top_srcdir)/src/misc
COMMON = ../../../Workspace.cpp ../../../Toolbox.cpp

TESTS = t_dft
check_PROGRAMS = $(TESTS)

t_dft_SOURCES = $(COMMON) t_dft.cpp
t_dft_LDFLAGS =  $(BLAS) $(LAPACK)

