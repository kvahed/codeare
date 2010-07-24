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
SERVER_SIDE = ReconContext.h ReconContext.cpp ReconStrategy.h ReconStrategy.cpp ReconServant.h ReconServant.cpp reconserver.cpp 

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE)

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE)
reconserver_LDFLAGS = -dlopen DummyRecon.la -dlopen DumpToFile.la -dlopen GenericRegrid.la -dlopen HannWindow.la -dlopen InvertOrder.la -dlopen MedianFilter.la -dlopen MedianFilter_OMP.la

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la HannWindow.la InvertOrder.la MedianFilter.la MedianFilter_OMP.la 

DummyRecon_la_SOURCES = DummyRecon.h DummyRecon.cpp
DummyRecon_la_LDFLAGS = -module

GenericRegrid_la_SOURCES = GenericRegrid.h GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = -module

MedianFilter_la_SOURCES = MedianFilter.h MedianFilter.cpp
MedianFilter_la_LDFLAGS = -module

MedianFilter_OMP_la_SOURCES = MedianFilter_OMP.h MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = -module

HannWindow_la_SOURCES = HannWindow.h HannWindow.cpp
HannWindow_la_LDFLAGS = -module

InvertOrder_la_SOURCES = InvertOrder.h InvertOrder.cpp
InvertOrder_la_LDFLAGS = -module

DumpToFile_la_SOURCES = DumpToFile.h DumpToFile.cpp
DumpToFile_la_LDFLAGS = -module

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc
