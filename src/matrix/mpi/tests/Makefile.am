AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/misc
COMMON = ../../../Workspace.cpp ../../../Toolbox.cpp

TESTS = t_alloc
check_PROGRAMS = $(TESTS)

t_alloc_SOURCES = $(COMMON) t_alloc.cpp
t_alloc_LDFLAGS = $(BLAS) $(LAPACK) $(SCALAPACK_LIBS) 
