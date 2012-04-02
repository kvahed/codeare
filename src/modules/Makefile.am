# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I../matrix -I../nfft -I../tinyxml -I../misc

EXTRA_DIST=GPUSimulator.cl

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la LapackTests.la SpatialDomain.la KTPoints.la RelativeSensitivities.la GRAPPA.la DirectMethod.la AFI.la

if USE_NFFT3
    lib_LTLIBRARIES += NuFFT.la CGSENSE.la NuFFT_OMP.la 
endif

if USE_GSL
    lib_LTLIBRARIES += CompressedSensing.la
endif

# Common linker flags
COMMON_la_LDFLAGS = -module -avoid-version @PLUGIN_LDFLAGS@

# Libraries
DummyRecon_la_SOURCES = DummyRecon.hpp DummyRecon.cpp
DummyRecon_la_LDFLAGS = $(COMMON_la_LDFLAGS)
DummyRecon_la_LIBADD = $(COMMON_la_LIBADD)

GenericRegrid_la_SOURCES = GenericRegrid.hpp GenericRegrid.cpp
GenericRegrid_la_LDFLAGS = $(COMMON_la_LDFLAGS)
GenericRegrid_la_LIBADD = $(COMMON_la_LIBADD)

MedianFilter_OMP_la_SOURCES = MedianFilter_OMP.hpp MedianFilter_OMP.cpp
MedianFilter_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS)
MedianFilter_OMP_la_LIBADD = $(COMMON_la_LIBADD)

NuFFT_la_SOURCES = NuFFT.hpp NuFFT.cpp
NuFFT_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3
NuFFT_la_LIBADD = $(COMMON_la_LIBADD)

NuFFT_OMP_la_SOURCES = NuFFT_OMP.hpp NuFFT_OMP.cpp OMP.hpp
NuFFT_OMP_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3
NuFFT_OMP_la_LIBADD = $(COMMON_la_LIBADD)

DumpToFile_la_SOURCES = DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS)
DumpToFile_la_LIBADD = $(COMMON_la_LIBADD)

CGSENSE_la_SOURCES = CGSENSE.hpp CGSENSE.cpp Noise.hpp ../matrix/NCSENSE.hpp ../matrix/SEM.hpp OMP.hpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 $(BLAS) $(LAPACK)
CGSENSE_la_LIBADD = $(COMMON_la_LIBADD)

SystemCmd_la_SOURCES = SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)
SystemCmd_la_LIBADD = $(COMMON_la_LIBADD)

LapackTests_la_SOURCES = LapackTests.hpp LapackTests.cpp
LapackTests_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
LapackTests_la_LIBADD = $(COMMON_la_LIBADD)

GRAPPA_la_SOURCES = GRAPPA.hpp GRAPPA.cpp
GRAPPA_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
GRAPPA_la_LIBADD = $(COMMON_la_LIBADD)

SpatialDomain_la_SOURCES = SpatialDomain.hpp SpatialDomain.cpp
SpatialDomain_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
SpatialDomain_la_LIBADD = $(COMMON_la_LIBADD)

KTPoints_la_SOURCES = KTPoints.hpp KTPoints.cpp
KTPoints_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
KTPoints_la_LIBADD = $(COMMON_la_LIBADD)

RelativeSensitivities_la_SOURCES = RelativeSensitivities.hpp RelativeSensitivities.cpp
RelativeSensitivities_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG)
RelativeSensitivities_la_LIBADD = $(COMMON_la_LIBADD) $(INSIGHT_LIB_FLAG) 
RelativeSensitivities_la_CXXFLAGS = $(AM_CXXFLAGS) $(INC_INSIGHT_H)

AFI_la_SOURCES = AFI.hpp AFI.cpp
AFI_la_LDFLAGS = $(COMMON_la_LDFLAGS)
AFI_la_LIBADD = $(COMMON_la_LIBADD)

CompressedSensing_la_SOURCES = CompressedSensing.hpp CompressedSensing.cpp
CompressedSensing_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
CompressedSensing_la_LIBADD = $(COMMON_la_LIBADD) 

DirectMethod_la_SOURCES = DirectMethod.hpp DirectMethod.cpp SimulationStrategy.hpp SimulationContext.hpp CPUSimulator.cpp CPUSimulator.hpp
if USE_OPENCL
    DirectMethod_la_SOURCES += GPUSimulator.cpp GPUSimulator.hpp
endif
DirectMethod_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
DirectMethod_la_LIBADD = $(COMMON_la_LIBADD)

CLEANFILES = RRSModule.hh RRSModuleSK.cc

