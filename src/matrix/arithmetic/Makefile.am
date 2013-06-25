AM_CXXFLAGS = -I../ -I../linalg -I../../ -I../../misc -I../io
COMMON = ../Matrix.hpp Trigonometry.hpp AVXTraits.hpp SSE2Traits.hpp SIMD.hpp SIMDTraits.hpp

TESTS  = t_sin t_cos t_tan t_exp t_asin t_acos t_atan t_atan2 t_sinh t_cosh t_tanh t_log t_flipud t_fliplr
TESTS += t_log10 t_emul t_ediv t_eadd t_esub t_nrmse t_epow t_floor t_round t_ceil t_size t_resize t_isnan 
TESTS += #t_gtoreq

check_PROGRAMS = $(TESTS)

t_sin_SOURCES = $(COMMON) t_sin.cpp
t_cos_SOURCES = $(COMMON) t_cos.cpp
t_tan_SOURCES = $(COMMON) t_tan.cpp
t_asin_SOURCES = $(COMMON) t_asin.cpp
t_acos_SOURCES = $(COMMON) t_acos.cpp
t_atan_SOURCES = $(COMMON) t_atan.cpp
t_atan2_SOURCES = $(COMMON) t_atan2.cpp
t_sinh_SOURCES = $(COMMON) t_sinh.cpp
t_cosh_SOURCES = $(COMMON) t_cosh.cpp
t_tanh_SOURCES = $(COMMON) t_tanh.cpp
t_exp_SOURCES = $(COMMON) t_exp.cpp
t_log_SOURCES = $(COMMON) t_log.cpp
t_log10_SOURCES = $(COMMON) t_log10.cpp
t_emul_SOURCES = $(COMMON) t_emul.cpp
t_ediv_SOURCES = $(COMMON) t_ediv.cpp
t_eadd_SOURCES = $(COMMON) t_eadd.cpp
t_esub_SOURCES = $(COMMON) t_esub.cpp
t_nrmse_SOURCES = $(COMMON) t_nrmse.cpp
t_epow_SOURCES = $(COMMON) t_epow.cpp
t_floor_SOURCES = $(COMMON) t_floor.cpp
t_round_SOURCES = $(COMMON) t_round.cpp
t_ceil_SOURCES = $(COMMON) t_ceil.cpp
t_size_SOURCES = $(COMMON) t_size.cpp
#t_gtoreq_SOURCES = $(COMMON) t_gtoreq.cpp
t_flipud_SOURCES = $(COMMON) t_flipud.cpp
t_fliplr_SOURCES = $(COMMON) t_fliplr.cpp
t_resize_SOURCES = $(COMMON) t_resize.cpp
t_isnan_SOURCES = $(COMMON) t_isnan.cpp


