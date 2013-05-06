AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/dwt -I$(top_srcdir)/src/misc
COMMON = ../../../Workspace.cpp ../../../Toolbox.cpp

TESTS = t_dwt
check_PROGRAMS = $(TESTS)

t_dwt_SOURCES = $(COMMON) t_dwt.cpp
t_dwt_LDFLAGS =  $(BLAS) $(LAPACK)

