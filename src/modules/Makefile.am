# Add idl to dependencies.
../%.hh: ../%.idl
	omniidl -bcxx $<
../%SK.cc: ../%.idl
	omniidl -bcxx $<

AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/linalg -I$(top_srcdir)/src/matrix/ft -I$(top_srcdir)/src/misc

EXTRA_DIST = GPUSimulator.cl

# Rules for recon strategies
lib_LTLIBRARIES = DummyRecon.la DumpToFile.la GenericRegrid.la MedianFilter_OMP.la SystemCmd.la SpatialDomain.la KTPoints.la RelativeSensitivities.la GRAPPA.la DirectMethod.la AFI.la SENSE.la

if USE_NFFT3
    lib_LTLIBRARIES += NuFFT.la CGSENSE.la 
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

DumpToFile_la_SOURCES = DumpToFile.hpp DumpToFile.cpp
DumpToFile_la_LDFLAGS = $(COMMON_la_LDFLAGS)
DumpToFile_la_LIBADD = $(COMMON_la_LIBADD)

SENSE_la_SOURCES = SENSE.hpp SENSE.cpp ../matrix/CSENSE.hpp
SENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
SENSE_la_LIBADD = $(COMMON_la_LIBADD)

CGSENSE_la_SOURCES = CGSENSE.hpp CGSENSE.cpp Noise.hpp ../matrix/NCSENSE.hpp ../matrix/SEM.hpp
CGSENSE_la_LDFLAGS = $(COMMON_la_LDFLAGS) -lnfft3 $(BLAS) $(LAPACK)
CGSENSE_la_LIBADD = $(COMMON_la_LIBADD)

SystemCmd_la_SOURCES = SystemCmd.hpp SystemCmd.cpp
SystemCmd_la_LDFLAGS = $(COMMON_la_LDFLAGS)
SystemCmd_la_LIBADD = $(COMMON_la_LIBADD)

GRAPPA_la_SOURCES = ../matrix/ft/GRAPPA.f ../matrix/ft/CGRAPPA.hpp GRAPPA.hpp GRAPPA.cpp
GRAPPA_la_DEPENDENCIES = 
GRAPPA_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
GRAPPA_la_LIBADD = $(COMMON_la_LIBADD)

SpatialDomain_la_SOURCES = SpatialDomain.hpp SpatialDomain.cpp
SpatialDomain_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
SpatialDomain_la_LIBADD = $(COMMON_la_LIBADD)

KTPoints_la_SOURCES = KTPoints.hpp KTPoints.cpp
KTPoints_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
KTPoints_la_LIBADD = $(COMMON_la_LIBADD)

RelativeSensitivities_la_SOURCES = RelativeSensitivities.hpp RelativeSensitivities.cpp
RelativeSensitivities_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK) $(INSIGHT_LIB_DIR_FLAG) -llapack
RelativeSensitivities_la_LIBADD = $(COMMON_la_LIBADD) $(INSIGHT_LIB_FLAG) 
RelativeSensitivities_la_CXXFLAGS = $(AM_CXXFLAGS) $(INC_INSIGHT_H)

AFI_la_SOURCES = AFI.hpp AFI.cpp
AFI_la_LDFLAGS = $(COMMON_la_LDFLAGS)
AFI_la_LIBADD = $(COMMON_la_LIBADD)

CompressedSensing_la_SOURCES = CompressedSensing.hpp CompressedSensing.cpp
CompressedSensing_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(top_srcdir)/src/matrix/dwt
CompressedSensing_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
CompressedSensing_la_LIBADD = $(COMMON_la_LIBADD) 

DirectMethod_la_SOURCES = DirectMethod.hpp DirectMethod.cpp SimulationStrategy.hpp SimulationContext.hpp CPUSimulator.cpp CPUSimulator.hpp #GPUSimulator.cpp GPUSimulator.hpp
DirectMethod_la_LDFLAGS = $(COMMON_la_LDFLAGS) $(BLAS) $(LAPACK)
DirectMethod_la_LIBADD = $(COMMON_la_LIBADD)

CLEANFILES = RRSModule.hh RRSModuleSK.cc

