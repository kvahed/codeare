# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<
COMMON = ../tinyxml.h ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Rules for recon strategies
lib_LTLIBRARIES = CGSENSE.la DummyRecon.la DumpToFile.la GenericRegrid.la InvertOrder.la MedianFilter_OMP.la NuFFT.la 

DummyRecon_la_SOURCES = $(COMMON) DummyRecon.h DummyRecon.cpp
DummyRecon_la_LDFLAGS = -module

GenericRegrid_la_SOURCES = $(COMMON) GenericRegrid.h GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = -module

MedianFilter_OMP_la_SOURCES = $(COMMON) MedianFilter_OMP.h MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = -module

NuFFT_la_SOURCES = $(COMMON) NuFFT.h NuFFT.cpp
NuFFT_la_LDFLAGS = -module -lnoncart

InvertOrder_la_SOURCES = $(COMMON) InvertOrder.h InvertOrder.cpp
InvertOrder_la_LDFLAGS = -module

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.h DumpToFile.cpp
DumpToFile_la_LDFLAGS = -module
