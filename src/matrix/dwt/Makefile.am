SUBDIRS = tests

AM_CPPFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/tinyxml

lib_LTLIBRARIES = libcodeare-dwt.la
libcodeare_dwt_la_SOURCES = DWT.hpp DWT.cpp
libcodeare_dwt_la_LDFLAGS = -module -avoid-version @PLUGIN_LDFLAGS@

include_HEADERS = codeare-dwt.h
