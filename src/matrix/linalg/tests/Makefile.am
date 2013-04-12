AM_CXXFLAGS = -DNVERBOSE -I../ -I../../ -I../../../ -I../../../misc
COMMON = ../Matrix.hpp ../Lapack.hpp

TESTS = t_chol t_svd t_inv t_pinv t_gemm t_gemv t_norm t_dot t_eig
check_PROGRAMS = $(TESTS)

t_svd_SOURCES = $(COMMON) t_svd.cpp
t_svd_LDFLAGS =  $(BLAS) $(LAPACK)
t_chol_SOURCES = $(COMMON) t_chol.cpp
t_chol_LDFLAGS =  $(BLAS) $(LAPACK)
t_dot_SOURCES = $(COMMON) t_dot.cpp
t_dot_LDFLAGS =  $(BLAS) $(LAPACK)
t_eig_SOURCES = $(COMMON) t_eig.cpp
t_eig_LDFLAGS =  $(BLAS) $(LAPACK)
t_gemm_SOURCES = $(COMMON) t_gemm.cpp
t_gemm_LDFLAGS =  $(BLAS) $(LAPACK)
t_gemv_SOURCES = $(COMVON) t_gemv.cpp
t_gemv_LDFLAGS =  $(BLAS) $(LAPACK)
t_inv_SOURCES = $(COMMON) t_inv.cpp
t_inv_LDFLAGS =  $(BLAS) $(LAPACK)
t_norm_SOURCES = $(COMMON) t_norm.cpp
t_norm_LDFLAGS =  $(BLAS) $(LAPACK)
t_pinv_SOURCES = $(COMMON) t_pinv.cpp
t_pinv_LDFLAGS =  $(BLAS) $(LAPACK)
