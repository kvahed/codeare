SUBDIRS = tests

AM_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix/io -I$(top_srcdir)/src/matrix/arithmetic
lib_LTLIBRARIES = libcodeare-linalg.la
libcodeare_linalg_la_SOURCES = Lapack.cpp

include_HEADERS = codeare-linalg.h
