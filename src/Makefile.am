# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = lclient
if NETWORKING
bin_PROGRAMS += rclient reconserver
endif

AM_CXXFLAGS = -DTIXML_USE_STL -Imatrix -Imatrix/linalg -Itinyxml -Imisc $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl Configurable.hpp DllExport.h matrix/Matrix.hpp options.h options.cpp Workspace.hpp Workspace.cpp FunctorContainer.hpp FunctorContainer.cpp ReconStrategy.hpp ReconContext.hpp ReconContext.cpp Loader.hpp Toolbox.hpp Toolbox.cpp Connector.hpp Grid.hpp Grid.cpp Params.hpp

# OmniORB interface
NETWORKING_COMMON = RRSModule.hh RRSModuleSK.cc 

# Rules for clients
noinst_LTLIBRARIES = liblocal.la
if NETWORKING
noinst_LTLIBRARIES += libremote.la 
endif

# Local client (lclient)
liblocal_la_SOURCES = LocalConnector.hpp LocalConnector.cpp testclt.cpp
liblocal_la_CXXFLAGS = $(AM_CXXFLAGS) -DLOCAL
lclient_SOURCES = $(COMMON)
lclient_LDFLAGS = liblocal.la tinyxml/libtinyxml.la $(BLAS) $(LAPACK)  $(SCALAPACK_LIBS) @PLUGIN_LDFLAGS@
# -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/KTPoints.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la -dlopen modules/DirectMethod.la 
lclient_DEPENDENCIES = liblocal.la 
if USE_NFFT3
    lclient_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la
endif

# Remote client (rclient)
libremote_la_SOURCES = RemoteConnector.hpp RemoteConnector.cpp testclt.cpp
libremote_la_CXXFLAGS = $(AM_CXXFLAGS) -DREMOTE
rclient_SOURCES = $(COMMON) $(NETWORKING_COMMON)
rclient_LDFLAGS = $(SCALAPACK_LIBS) libremote.la tinyxml/libtinyxml.la $(BLAS) $(LAPACK) @PLUGIN_LDFLAGS@ -lomniORB4 -lomnithread -lomniDynamic4 
rclient_DEPENDENCIES = libremote.la

# Rules for reconserver
reconserver_SOURCES = $(COMMON) $(NETWORKING_COMMON) ReconServant.hpp ReconServant.cpp reconserver.cpp
reconserver_LDFLAGS = $(SCALAPACK_LIBS) @PLUGIN_LDFLAGS@  -lomniORB4 -lomnithread -lomniDynamic4 tinyxml/libtinyxml.la
#-dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/KTPoints.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la -dlopen modules/DirectMethod.la 
if USE_NFFT3
    reconserver_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la
endif

EXTRA_DIST = RRSModule.idl Makefile.win32 common.h
CLEANFILES = RRSModule.hh RRSModuleSK.cc

SUBDIRS = tinyxml modules matrix 
