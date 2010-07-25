# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = testclt reconserver

# Common sources
COMMON =  common.h Matrix.h Matrix.cpp options.h options.cpp RRSModule.hh RRSModuleSK.cc

# Client side stuff
CLIENT_SIDE =  ReconClient.h ReconClient.cpp testclt.cpp

# Reconstruction strategies & server side stuff
SERVER_SIDE = ReconContext.h ReconContext.cpp ReconStrategy.h ReconStrategy.cpp ReconServant.h ReconServant.cpp reconserver.cpp Loader.h

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE)

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE)
reconserver_LDFLAGS = -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/HannWindow.la -dlopen modules/InvertOrder.la -dlopen modules/MedianFilter.la -dlopen modules/MedianFilter_OMP.la

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc

SUBDIRS = modules
