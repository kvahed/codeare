# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = testclt reconserver

# Common sources
COMMON =  common.h Matrix.h Matrix.cpp options.h options.cpp RRSModule.hh RRSModuleSK.cc

# Client side stuff
CLIENT_SIDE =  ReconClient.h ReconClient.cpp testclt.cpp

# Reconstruction strategies & server side stuff
RECON_STRATEGIES = DummyRecon.h DummyRecon.cpp DumpToFile.h DumpToFile.cpp GenericRegrid.h GenericRegrid.cpp HannWindow.h HannWindow.cpp InvertOrder.h InvertOrder.cpp MedianFilter.h MedianFilter.cpp MedianFilter_OMP.h MedianFilter_OMP.cpp
SERVER_SIDE = ReconContext.h ReconContext.cpp ReconStrategy.h ReconStrategy.cpp ReconServant.h ReconServant.cpp reconserver.cpp 

# Rules for the binaries
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE)
reconserver_SOURCES = $(COMMON) $(RECON_STRATEGIES) $(SERVER_SIDE)

EXTRA_DIST = RRSModule.idl
