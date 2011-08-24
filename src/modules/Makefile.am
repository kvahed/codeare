# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS)

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la NuFFT_OMP.la NuFFT.la CGSENSE.la SpatialDomain.la RelativeSensitivities.la GRAPPA.la 

COMMON = ../tinyxml/tinystr.cpp ../tinyxml/tinystr.h ../tinyxml/tinyxml.cpp ../tinyxml/tinyxml.h ../tinyxml/tinyxmlerror.cpp ../tinyxml/tinyxmlparser.cpp
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

NuFFT_OMP_la_SOURCES = $(COMMON) NuFFT_OMP.hpp NuFFT_OMP.cpp nfftstub.h nfftstub.c OMP.hpp
NuFFT_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart

DumpToFile_la_SOURCES = $(COMMON) DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS) #$(HDF5)

CGSENSE_la_SOURCES = $(COMMON) CGSENSE.hpp CGSENSE.cpp nfftstub.h nfftstub.c Noise.hpp SEM.hpp OMP.hpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnoncart $(BLAS) $(LAPACK) #-Wl,-rpath,/Applications/MATLAB_R2011a.app/bin/maci64

SystemCmd_la_SOURCES = $(COMMON) SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)

LapackTests_la_SOURCES = $(COMMON) LapackTests.hpp LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

GRAPPA_la_SOURCES = $(COMMON) GRAPPA.hpp GRAPPA.cpp
GRAPPA_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(PYTHON_LDFLAGS) -lpython

SpatialDomain_la_SOURCES = $(COMMON) SpatialDomain.hpp SpatialDomain.cpp
SpatialDomain_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) 

RelativeSensitivities_la_SOURCES = $(COMMON) RelativeSensitivities.hpp RelativeSensitivities.cpp
RelativeSensitivities_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) 

CLEANFILES = RRSModule.hh RRSModuleSK.cc
