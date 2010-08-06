# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<
COMMON = #../RRSModule.hh ../RRSModuleSK.cc

AM_CPPFLAGS = -DTIXML_USE_STL

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la HannWindow.la InvertOrder.la MedianFilter_OMP.la 

DummyRecon_la_SOURCES = $(COMMON) DummyRecon.h DummyRecon.cpp
DummyRecon_la_LDFLAGS = -module

GenericRegrid_la_SOURCES = $(COMMON) GenericRegrid.h GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = -module

MedianFilter_OMP_la_SOURCES = $(COMMON) MedianFilter_OMP.h MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = -module

HannWindow_la_SOURCES = $(COMMON) HannWindow.h HannWindow.cpp
HannWindow_la_LDFLAGS = -module

InvertOrder_la_SOURCES = $(COMMON) InvertOrder.h InvertOrder.cpp
InvertOrder_la_LDFLAGS = -module

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.h DumpToFile.cpp
DumpToFile_la_LDFLAGS = -module
