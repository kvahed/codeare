# Add idl to dependencies.
%.hh: common.h %.idl
	omniidl -bcxx $<

%SK.cc: common.h %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = testclt reconserver

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl common.h Configurable.h DllExport.h Matrix.h options.h options.cpp RRSModule.hh RRSModuleSK.cc tinyxml.h tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp 

# Client side stuff
CLIENT_SIDE =  ReconClient.h ReconClient.cpp

# Reconstruction strategies & server side stuff
SERVER_SIDE =  Loader.h ReconContext.h ReconContext.cpp ReconServant.h ReconServant.cpp ReconStrategy.h ReconStrategy.cpp

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE) testclt.cpp

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE) reconserver.cpp
reconserver_LDFLAGS = -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la
if NONCART
  reconserver_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la
endif

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc

SUBDIRS = modules
