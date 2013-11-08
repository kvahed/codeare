AM_CXXFLAGS = -DTIXML_USE_STL $(OPENMP_CXXFLAGS) -I$(top_srcdir)/src/misc -I$(top_srcdir)/src/matrix -I$(top_srcdir)/src/matrix/io -I$(top_srcdir)/src/matrix/interp

TESTS = t_medfilt2 t_ppval t_interp1 t_resample
check_PROGRAMS = $(TESTS)

t_medfilt2_SOURCES = $(COMMON) t_medfilt2.cpp
t_ppval_SOURCES= $(COMMON) t_ppval.cpp
t_interp1_SOURCES= $(COMMON) t_interp1.cpp
t_resample_SOURCES= $(COMMON) t_resample.cpp
t_resample_CXXFLAGS = $(AM_CXXFLAGS) $(INC_INSIGHT_H)
t_resample_LDADD = $(COMMON_la_LDADD) $(INSIGHT_LIB_FLAG) $(INSIGHT_LIB_DIR_FLAG)



