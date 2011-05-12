# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

COMMON = ../tinyxml.h ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp ../Toolbox.hpp ../Toolbox.cpp

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la CGSENSE.la NuFFT.la

# Common linker flags
COMMON_la_LDFLAGS = -module

# Libraries
DummyRecon_la_SOURCES = $(COMMON) DummyRecon.hpp DummyRecon.cpp
DummyRecon_la_LDFLAGS = $(COMMON_la_LDFLAGS)

GenericRegrid_la_SOURCES = $(COMMON) GenericRegrid.hpp GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = $(COMMON_la_LDFLAGS)

MedianFilter_OMP_la_SOURCES = $(COMMON) MedianFilter_OMP.hpp MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS)

NuFFT_la_SOURCES = $(COMMON) NuFFT.hpp NuFFT.cpp nfftstub.h nfftstub.c
NuFFT_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(HDF5)

CGSENSE_la_SOURCES = $(COMMON) CGSENSE.hpp CGSENSE.cpp nfftstub.h nfftstub.c
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart

SystemCmd_la_SOURCES = $(COMMON) SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)

LapackTests_la_SOURCES = $(COMMON) LapackTests.hpp LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

CLEANFILES = RRSModule.hh RRSModuleSK.cc
