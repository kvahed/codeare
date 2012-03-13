# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = reconserver testclt 

AM_CXXFLAGS = -DTIXML_USE_STL -Imatrix -Itinyxml -Imisc $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl Configurable.h DllExport.h matrix/Matrix.hpp options.h options.cpp RRSModule.hh RRSModuleSK.cc DataBase.hpp DataBase.cpp FunctorContainer.hpp FunctorContainer.cpp ReconStrategy.hpp ReconContext.hpp ReconContext.cpp Loader.hpp Toolbox.hpp Toolbox.cpp

# Client side stuff
CLIENT_SIDE = RemoteConnector.hpp RemoteConnector.cpp Connector.hpp LocalConnector.hpp LocalConnector.cpp

# Reconstruction strategies & server side stuff
SERVER_SIDE = ReconServant.hpp ReconServant.cpp

# Rules for testclt
testclt_SOURCES = $(COMMON) $(CLIENT_SIDE)  testclt.cpp 
testclt_LDFLAGS = tinyxml/libtinyxml.la matrix/libmatrix.la $(BLAS) $(LAPACK) @PLUGIN_LDFLAGS@ 
testclt_DEPENDENCIES = matrix/libmatrix.la
#testclt_LDADD = 

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(SERVER_SIDE) reconserver.cpp
reconserver_LDFLAGS = $(HDF5) @PLUGIN_LDFLAGS@ -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/KTPoints.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la -dlopen modules/DirectMethod.la  tinyxml/libtinyxml.la matrix/libmatrix.la
reconserver_DEPENDENCIES = matrix/libmatrix.la
if USE_NFFT3
    reconserver_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la -dlopen modules/NuFFT_OMP.la 
endif

EXTRA_DIST = RRSModule.idl Makefile.win32
CLEANFILES = RRSModule.hh RRSModuleSK.cc common.hh commonSK.cc

SUBDIRS = tinyxml matrix nfft modules
