# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = testclt reconserver

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl common.h Configurable.h DllExport.h Matrix.h options.h options.cpp RRSModule.hh RRSModuleSK.cc tinyxml.h tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp Toolbox.cpp

# Client side stuff
CLIENT_SIDE =  ReconClient.h ReconClient.cpp ReconContext.h ReconContext.cpp  ReconStrategy.hpp

# Reconstruction strategies & server side stuff
SERVER_SIDE =  Loader.h ReconContext.h ReconContext.cpp ReconServant.h ReconServant.cpp ReconStrategy.hpp

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE) testclt.cpp
testclt_LDFLAGS = $(BLAS) $(LAPACK) 

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE) reconserver.cpp
reconserver_LDFLAGS = $(HDF5) -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc common.hh commonSK.cc

SUBDIRS = modules
