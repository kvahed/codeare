# Add idl to dependencies.
%.hh: %.idl
	omniidl -bcxx $<

%SK.cc: %.idl
	omniidl -bcxx $<

# Binaries 
bin_PROGRAMS = lclient
if NETWORKING
bin_PROGRAMS += rclient codeared
endif

AM_CXXFLAGS = -DTIXML_USE_STL -Icore -Imatrix -Imatrix/simd -Imatrix/linalg -Itinyxml -Imisc -Imatrix/arithmetic -Imatrix/interp -Imatrix/io -Imodules $(OPENMP_CXXFLAGS)

# Common sources
COMMON =  RRSModule.idl Configurable.hpp DllExport.h options.h options.cpp ReconStrategy.hpp Loader.hpp Connector.hpp Grid.hpp Grid.cpp

COMMON_LD_FLAGS=tinyxml/libtinyxml.la core/libcore.la

# OmniORB interface
NETWORKING_COMMON = RRSModule.hh RRSModuleSK.cc 

# Rules for clients
noinst_LTLIBRARIES = liblocal.la
if NETWORKING
noinst_LTLIBRARIES += libremote.la 
endif

# Local client (lclient)
liblocal_la_SOURCES = Loader.hpp LocalConnector.hpp LocalConnector.cpp testclt.cpp
liblocal_la_CXXFLAGS = $(AM_CXXFLAGS) -DLOCAL
lclient_SOURCES = $(COMMON)
lclient_LDFLAGS = $(COMMON_LD_FLAGS) liblocal.la $(BLAS) $(LAPACK)  $(SCALAPACK_LIBS) @PLUGIN_LDFLAGS@
# -dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/KTPoints.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la -dlopen modules/DirectMethod.la 
lclient_DEPENDENCIES = liblocal.la
if USE_NFFT3
    lclient_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la
endif

# Remote client (rclient)
libremote_la_SOURCES = Loader.hpp RemoteConnector.hpp RemoteConnector.cpp testclt.cpp
libremote_la_CXXFLAGS = $(AM_CXXFLAGS) -DREMOTE
rclient_SOURCES = $(COMMON) $(NETWORKING_COMMON)
rclient_LDFLAGS = $(SCALAPACK_LIBS) libremote.la $(COMMON_LD_FLAGS) $(BLAS) $(LAPACK) @PLUGIN_LDFLAGS@ -lomniORB4 -lomnithread -lomniDynamic4 
rclient_DEPENDENCIES = libremote.la

# Rules for codeared
codeared_SOURCES = Loader.hpp $(COMMON) $(NETWORKING_COMMON) ReconServant.hpp ReconServant.cpp codeared.cpp
codeared_CXXFLAGS = $(AM_CXXFLAGS) -Imongoose
codeared_LDFLAGS = $(SCALAPACK_LIBS) @PLUGIN_LDFLAGS@  -lomniORB4 -lomnithread -lomniDynamic4 $(COMMON_LD_FLAGS) mongoose/libmongoose.la
#-dlopen modules/DummyRecon.la -dlopen modules/DumpToFile.la -dlopen modules/GenericRegrid.la -dlopen modules/SystemCmd.la -dlopen modules/MedianFilter_OMP.la -dlopen modules/SpatialDomain.la -dlopen modules/KTPoints.la -dlopen modules/RelativeSensitivities.la -dlopen modules/GRAPPA.la -dlopen modules/DirectMethod.la 
if USE_NFFT3
    codeared_LDFLAGS += -dlopen modules/CGSENSE.la -dlopen modules/NuFFT.la
endif

EXTRA_DIST = RRSModule.idl Makefile.win32 common.h
CLEANFILES = RRSModule.hh RRSModuleSK.cc

SUBDIRS = tinyxml core modules matrix mongoose
