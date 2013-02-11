AM_CXXFLAGS = -DVERBOSE -I../ -I../../ -I../../misc
COMMON = ../Matrix.hpp ../Lapack.hpp

TESTS = t_svd t_eig 
check_PROGRAMS = $(TESTS)

t_svd_SOURCES = $(COMMON) t_svd.cpp
t_svd_LDFLAGS =  $(BLAS) $(LAPACK)
t_eig_SOURCES = $(COMMON) t_eig.cpp
t_eig_LDFLAGS =  $(BLAS) $(LAPACK)
#t_cos_SOURCES = $(COMMON) t_cos.cpp
#t_tan_SOURCES = $(COMMON) t_tan.cpp
#t_asin_SOURCES = $(COMMON) t_asin.cpp
#t_acos_SOURCES = $(COMMON) t_acos.cpp
#t_atan_SOURCES = $(COMMON) t_atan.cpp
#t_atan2_SOURCES = $(COMMON) t_atan2.cpp
#t_sinh_SOURCES = $(COMMON) t_sinh.cpp
#t_cosh_SOURCES = $(COMMON) t_cosh.cpp
#t_tanh_SOURCES = $(COMMON) t_tanh.cpp
#t_exp_SOURCES = $(COMMON) t_exp.cpp
#t_log_SOURCES = $(COMMON) t_log.cpp
#t_log10_SOURCES = $(COMMON) t_log10.cpp
