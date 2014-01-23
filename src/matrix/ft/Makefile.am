SUBDIRS = tests

AM_CPPFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/simd -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/io -I$(top_srcdir)/src/matrix/arithmetic -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/tinyxml

lib_LTLIBRARIES = libcodeare-ft.la
libcodeare_ft_la_SOURCES = DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp CGRAPPA.cpp CGRAPPA.hpp
libcodeare_ft_la_LDFLAGS = -module -avoid-version @PLUGIN_LDFLAGS@
if USE_NFFT3
    libcodeare_ft_la_SOURCES += NFFT.hpp NFFT.cpp NCSENSE.hpp NCSENSE.cpp
	libcodeare_ft_la_LIBADD += -lnfft3
endif

include_HEADERS = codeare-ft.h
