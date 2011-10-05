# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) $(INC_INSIGHT_H)

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la SpatialDomain.la RelativeSensitivities.la GRAPPA.la DirectMethod.la

if USE_NFFT3
    lib_LTLIBRARIES += NuFFT_OMP.la NuFFT.la CGSENSE.la 
endif

if USE_GSL
    lib_LTLIBRARIES += CompressedSensing.la
endif

COMMON = ../Toolbox.hpp ../Toolbox.cpp ../tinyxml/tinystr.cpp ../tinyxml/tinystr.h ../tinyxml/tinyxml.cpp ../tinyxml/tinyxml.h ../tinyxml/tinyxmlerror.cpp ../tinyxml/tinyxmlparser.cpp ../tinyxml/xpath_processor.h ../tinyxml/xpath_processor.cpp ../tinyxml/xpath_expression.cpp ../tinyxml/xpath_expression.h ../tinyxml/xpath_static.cpp ../tinyxml/xpath_static.h ../tinyxml/action_store.cpp ../tinyxml/action_store.h ../tinyxml/xpath_stack.cpp ../tinyxml/xpath_stack.h ../tinyxml/lex_util.cpp ../tinyxml/lex_util.h ../tinyxml/xpath_stream.cpp ../tinyxml/xpath_stream.h ../tinyxml/node_set.cpp ../tinyxml/node_set.h ../tinyxml/xml_util.cpp ../tinyxml/xml_util.h ../tinyxml/xpath_syntax.cpp ../tinyxml/xpath_syntax.h ../tinyxml/tokenlist.cpp ../tinyxml/tokenlist.h FFT.hpp FFT.cpp

# Common linker flags
COMMON_la_LDFLAGS = -module

# Libraries
DummyRecon_la_SOURCES = $(COMMON) DummyRecon.hpp DummyRecon.cpp
DummyRecon_la_LDFLAGS = $(COMMON_la_LDFLAGS)

GenericRegrid_la_SOURCES = $(COMMON) GenericRegrid.hpp GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = $(COMMON_la_LDFLAGS)

MedianFilter_OMP_la_SOURCES = $(COMMON) MedianFilter_OMP.hpp MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS)

NuFFT_la_SOURCES = $(COMMON) NuFFT.hpp NuFFT.cpp nfftstub.h nfftstub.c
NuFFT_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3

NuFFT_OMP_la_SOURCES = $(COMMON) NuFFT_OMP.hpp NuFFT_OMP.cpp nfftstub.h nfftstub.c OMP.hpp
NuFFT_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS)

CGSENSE_la_SOURCES = $(COMMON) CGSENSE.hpp CGSENSE.cpp nfftstub.h nfftstub.c Noise.hpp SEM.hpp OMP.hpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 $(BLAS) $(LAPACK)

SystemCmd_la_SOURCES = $(COMMON) SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)

LapackTests_la_SOURCES = $(COMMON) LapackTests.hpp LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

GRAPPA_la_SOURCES = $(COMMON) GRAPPA.hpp GRAPPA.cpp
GRAPPA_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

SpatialDomain_la_SOURCES = $(COMMON) SpatialDomain.hpp SpatialDomain.cpp
SpatialDomain_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

RelativeSensitivities_la_SOURCES = $(COMMON) RelativeSensitivities.hpp RelativeSensitivities.cpp
RelativeSensitivities_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG) $(INSIGHT_LIB_FLAG) 

CompressedSensing_la_SOURCES = $(COMMON) DWT.cpp DWT.hpp CompressedSensing.hpp CompressedSensing.cpp
CompressedSensing_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG) $(INSIGHT_LIB_FLAG) 

DirectMethod_la_SOURCES = $(COMMON) DirectMethod.hpp DirectMethod.cpp BlochSim.hpp
DirectMethod_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

CLEANFILES = RRSModule.hh RRSModuleSK.cc
