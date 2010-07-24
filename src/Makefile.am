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
#RECON_STRATEGIES = DummyRecon.h DummyRecon.cpp DumpToFile.h DumpToFile.cpp GenericRegrid.h GenericRegrid.cpp HannWindow.h HannWindow.cpp InvertOrder.h InvertOrder.cpp MedianFilter.h MedianFilter.cpp MedianFilter_OMP.h MedianFilter_OMP.cpp

lib_LTLIBRARIES = libDumpToFile.la libDummyRecon.la libGenericRegrid.la libHannWindow.la libInvertOrder.la libMedianFilter.la libMedianFilter_OMP.la  

libDumpToFile_la_SOURCES = DumpToFile.h DumpToFile.cpp
libDumpToFile_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libDummyRecon_la_SOURCES = DummyRecon.h DummyRecon.cpp
libDummyRecon_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libGenericRegrid_la_SOURCES = GenericRegrid.h GenericRegrid.cpp
libGenericRegrid_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libHannWindow_la_SOURCES = HannWindow.h HannWindow.cpp
libHannWindow_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libInvertOrder_la_SOURCES = InvertOrder.h InvertOrder.cpp
libInvertOrder_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libMedianFilter_la_SOURCES = MedianFilter.h MedianFilter.cpp
libMedianFilter_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

libMedianFilter_OMP_la_SOURCES = MedianFilter_OMP.h MedianFilter_OMP.cpp
libMedianFilter_OMP_la_LDFLAGS = $(all_libraries) --version-info 1.0.0 

SERVER_SIDE = ReconContext.h ReconStrategy.h ReconServant.h ReconServant.cpp reconserver.cpp 

# Rules for the binaries
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE)
#reconserver_SOURCES = $(COMMON) $(RECON_STRATEGIES) $(SERVER_SIDE)
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE)

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc
