#include <complex>

#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#ifdef __SSSE3__
#include <tmmintrin.h>
#endif 
#ifdef __SSSE41__
#include <smmintrin.h>
#endif 
#ifdef __SSSE42__
#include <nmmintrin.h>
#endif 
#ifdef __AVX__
#include <immintrin.h>
#endif

namespace SSE {

	template<class T> 
	struct SSETraits {};

#if defined __AVX__
#include "AVXTraits.hpp"
#elif defined __SSE2__
#include "SSE2Traits.hpp"
#endif

} // namespace SSE 
