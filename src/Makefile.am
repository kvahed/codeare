# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = testclt reconserver

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl common.h Configurable.h DllExport.h Matrix.h options.h options.cpp RRSModule.hh RRSModuleSK.cc Toolbox.cpp tinyxml/tinystr.cpp tinyxml/tinystr.h tinyxml/tinyxml.cpp tinyxml/tinyxml.h tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp tinyxml/xpath_processor.h tinyxml/xpath_processor.cpp tinyxml/xpath_expression.cpp tinyxml/xpath_expression.h tinyxml/xpath_static.cpp tinyxml/xpath_static.h tinyxml/action_store.cpp tinyxml/action_store.h tinyxml/xpath_stack.cpp tinyxml/xpath_stack.h tinyxml/lex_util.cpp tinyxml/lex_util.h tinyxml/xpath_stream.cpp tinyxml/xpath_stream.h tinyxml/node_set.cpp tinyxml/node_set.h tinyxml/xml_util.cpp tinyxml/xml_util.h tinyxml/xpath_syntax.cpp tinyxml/xpath_syntax.h tinyxml/tokenlist.cpp tinyxml/tokenlist.h

# Client side stuff
CLIENT_SIDE =  ReconClient.h ReconClient.cpp ReconContext.h ReconContext.cpp  ReconStrategy.hpp

# Reconstruction strategies & server side stuff
SERVER_SIDE =  Loader.h ReconContext.h ReconContext.cpp ReconServant.h ReconServant.cpp ReconStrategy.hpp

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE) testclt.cpp
testclt_LDFLAGS = $(BLAS) $(LAPACK) 

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE) reconserver.cpp
reconserver_LDFLAGS = $(HDF5) -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la
if USE_NFFT3
    reconserver_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la -dlopen modules/NuFFT_OMP.la 
endif

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc common.hh commonSK.cc

SUBDIRS = modules tinyxml
