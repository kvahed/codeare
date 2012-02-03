# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) $(INC_INSIGHT_H) -I../matrix -I../nfft

EXTRA_DIST=GPUSimulator.cl

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la SpatialDomain.la KTPoints.la RelativeSensitivities.la GRAPPA.la DirectMethod.la AFI.la

if USE_NFFT3
    lib_LTLIBRARIES += NuFFT_OMP.la NuFFT.la CGSENSE.la 
endif

if USE_GSL
    lib_LTLIBRARIES += CompressedSensing.la
endif

# Common linker flags
COMMON_la_LDFLAGS = -module -avoid-version @PLUGIN_LDFLAGS@ ../tinyxml/libtinyxml.a ../matrix/libmatrix.a

# Libraries
DummyRecon_la_SOURCES = DummyRecon.hpp DummyRecon.cpp
DummyRecon_la_LDFLAGS = $(COMMON_la_LDFLAGS)

GenericRegrid_la_SOURCES = GenericRegrid.hpp GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = $(COMMON_la_LDFLAGS)

MedianFilter_OMP_la_SOURCES = MedianFilter_OMP.hpp MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS)

NuFFT_la_SOURCES = NuFFT.hpp NuFFT.cpp
NuFFT_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 ../nfft/libnfftstub.a

NuFFT_OMP_la_SOURCES = NuFFT_OMP.hpp NuFFT_OMP.cpp OMP.hpp
NuFFT_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 ../nfft/libnfftstub.a

DumpToFile_la_SOURCES = DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS)

CGSENSE_la_SOURCES = CGSENSE.hpp CGSENSE.cpp Noise.hpp SEM.hpp OMP.hpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 $(BLAS) $(LAPACK) ../nfft/libnfftstub.a

SystemCmd_la_SOURCES = SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)

LapackTests_la_SOURCES = LapackTests.hpp LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

GRAPPA_la_SOURCES = GRAPPA.hpp GRAPPA.cpp
GRAPPA_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

SpatialDomain_la_SOURCES = SpatialDomain.hpp SpatialDomain.cpp
SpatialDomain_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

KTPoints_la_SOURCES = KTPoints.hpp KTPoints.cpp
KTPoints_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

RelativeSensitivities_la_SOURCES = RelativeSensitivities.hpp RelativeSensitivities.cpp
RelativeSensitivities_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG) $(INSIGHT_LIB_FLAG) 

AFI_la_SOURCES = AFI.hpp AFI.cpp
AFI_la_LDFLAGS = $(COMMON_la_LDFLAGS)

CompressedSensing_la_SOURCES = CompressedSensing.hpp CompressedSensing.cpp
CompressedSensing_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG) $(INSIGHT_LIB_FLAG) 

DirectMethod_la_SOURCES = DirectMethod.hpp DirectMethod.cpp SimulationStrategy.hpp SimulationContext.hpp CPUSimulator.cpp CPUSimulator.hpp GPUSimulator.cpp GPUSimulator.hpp 
DirectMethod_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)

CLEANFILES = RRSModule.hh RRSModuleSK.cc
