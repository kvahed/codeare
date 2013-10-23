AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/unwarp -I$(top_srcdir)/src/matrix/arithmetic -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix/io
COMMON = ../../../Workspace.cpp ../../../Toolbox.cpp

TESTS = t_unwrap
check_PROGRAMS = $(TESTS)

t_unwrap_SOURCES = $(COMMON) t_unwrap.cpp
t_unwrap_LDFLAGS =  $(BLAS) $(LAPACK)

