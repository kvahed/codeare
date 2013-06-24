SUBDIRS = io arithmetic linalg ft mpi
EXTRA_DIST = Access.hpp Algos.hpp CGRAPPA.hpp CSENSE.hpp CX.hpp Complex.hpp Creators.hpp DFT.hpp DWT.hpp DWTTraits.hpp FFTWTraits.hpp FT.hpp IO.hpp Interpolate.hpp Lapack.hpp LapackTraits.hpp MCGLS.hpp MITK.hpp MRI.hpp Math.hpp Matrix.hpp NCSENSE.hpp NFFT.hpp NFFTTraits.hpp OMP.hpp PIO.hpp PMatrix.hpp PolyVal.hpp SEM.hpp SSE.hpp SSETraits.hpp ScalapackTraits.hpp Statistics.hpp TVOP.hpp XProtocolParser.hpp 

AM_CPPFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/io  -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/tinyxml

lib_LTLIBRARIES = libcodeare-matrix.la
libcodeare_matrix_la_SOURCES = Matrix.hpp Matrix.cpp Algos.hpp
libcodeare_matrix_la_LDFLAGS = -module -avoid-version @PLUGIN_LDFLAGS@

include_HEADERS = codeare-matrix.h
