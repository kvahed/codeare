# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

COMMON = ../tinyxml.h ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la

if NONCART
  lib_LTLIBRARIES += CGSENSE.la NuFFT.la
endif

# Common linker flags
COMMON_la_LDFLAGS = -module $(BLAS) $(LAPACK)

# Libraries
DummyRecon_la_SOURCES = $(COMMON) DummyRecon.h DummyRecon.cpp
DummyRecon_la_LDFLAGS = $(COMMON_la_LDFLAGS)

GenericRegrid_la_SOURCES = $(COMMON) GenericRegrid.h GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = $(COMMON_la_LDFLAGS)

MedianFilter_OMP_la_SOURCES = $(COMMON) MedianFilter_OMP.h MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS)

NuFFT_la_SOURCES = $(COMMON) NuFFT.h NuFFT.cpp
NuFFT_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.h DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS)

CGSENSE_la_SOURCES = $(COMMON) CGSENSE.h CGSENSE.cpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart

SystemCmd_la_SOURCES = $(COMMON) SystemCmd.h SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)

LapackTests_la_SOURCES = $(COMMON) LapackTests.h LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS)

CLEANFILES = RRSModule.hh RRSModuleSK.cc
