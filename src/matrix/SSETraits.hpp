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

#ifdef __AVX__

	template<>
	struct SSETraits<double> {
		
		typedef __m256d Register;         /**< @brief register type */
		static const unsigned int ne = 4; /**< @brief # of processed elements */
		
		/**
		 * @brief     AVX load packed aligned
		 */
		static inline Register 
		loada (const double* p) {
			return _mm256_load_pd (p); 
		}
		
		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const double* p) {
			return _mm256_loadu_pd (p); 
		}
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (double* p, Register a) {
			_mm256_store_pd (p, a); 
		}
		
		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (double* p, Register a) {
			_mm256_storeu_pd (p, a); 
		}
		
		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm256_add_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm256_sub_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register &a, const Register &b) {
			return _mm256_mul_pd(a, b);
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		
		static inline Register 
		muls (const Register &a, const Register &b) {
			return a;//_mm256_mul_sd(a, b);
		}


		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &a, const Register &b) {
			return _mm256_div_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		/*
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm256_div_sd(a, b);
		}
		*/
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm256_sqrt_pd(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		/*
		static inline Register 
		    sqrts (const Register &a, const Register &b) {
		    return _mm256_sqrt_sd(a,b);
		}
		*/
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm256_min_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		/*
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm256_min_sd(a, b);
		}
		*/

		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm256_max_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		/*
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm256_max_sd(a, b);
		}
		*/
	}; // SSETraits<double>
	


	template<>
	struct SSETraits<float> {
		
		typedef __m256 Register;         /**< @brief register type */
		static const unsigned int ne = 8; /**< @brief # of processed elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const float* p) {
			return _mm256_load_ps (p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const float* p) {
			return _mm256_loadu_ps (p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (float* p, Register a) {
			_mm256_store_ps (p, a); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (float* p, Register a) {
			_mm256_storeu_ps (p, a); 
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm256_add_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */

		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm256_sub_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register &a, const Register &b) {
			return _mm256_mul_ps(a, b);
		}
	
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return a;//_mm256_mul_ps(a, b);
		}
	
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &a, const Register &b) {
			return _mm256_div_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		/*
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm256_div_ss(a, b);
		}
		*/
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm256_sqrt_ps(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		/*
		static inline Register 
		sqrts (const Register &a) {
			return _mm256_sqrt_ss(a);
		}
		*/
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm256_min_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		/*
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm256_min_ss(a, b);
		}
		*/
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm256_max_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		/*
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm256_max_ss(a, b);
		}
		*/
	}; // SSETraits<float>	

#else
	
	template<>
	struct SSETraits< std::complex<double> > {
		
		typedef __m128d Register;         /**< @brief register type */
		static const unsigned int ne = 1; /**< @brief # of processed elements */
		static const unsigned int ns = 2; /**< @brief # of sub elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const std::complex<double>* p) {
			return _mm_load_pd ((double*)p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const std::complex<double>* p) {
			return _mm_loadu_pd ((double*)p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loadoa (const std::complex<double>* p) {
			return _mm_load_pd ((double*)p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadou (const std::complex<double>* p) {
			return _mm_loadu_pd ((double*)p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (std::complex<double>* p, Register a) {
			_mm_store_pd ((double*)p, a); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (std::complex<double>* p, Register a) {
			_mm_storeu_pd ((double*)p, a); 
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm_add_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register &a, const Register &b) {
			return _mm_add_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm_sub_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register &a, const Register &b) {
			return _mm_sub_sd(a, b);
		}
		
		/**
		 * @brief     SSE3 packed multiplication
		 */
		static inline Register 
		mulp (const Register &a, const Register &b) {

			Register c, d, e, f;

			c = _mm_mul_pd (a,b);				
			d = _mm_shuffle_pd (c,c,0x1);		
			e = _mm_shuffle_pd (b,b,0x1);		
			f = _mm_sub_pd (c,d);				
			d = _mm_mul_pd (a,e);				
			e = _mm_shuffle_pd (d,d,0x1);		
			e = _mm_add_pd (d,e);	
			
			return _mm_shuffle_pd (f,e,0x2);

		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return mulp(a, b);
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &a, const Register &b) {
			return _mm_div_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm_div_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm_sqrt_pd(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register &a, const Register &b) {
			return _mm_sqrt_sd(a,b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm_min_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm_min_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm_max_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm_max_sd(a, b);
		}
		
	}; // SSETraits< std::complex<double> >

	template<>
	struct SSETraits< double > {
		
		typedef __m128d Register;         /**< @brief register type */
		static const unsigned int ne = 2; /**< @brief # of processed elements */
		static const unsigned int ns = 1; /**< @brief # of sub elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const double* p) {
			return _mm_load_pd (p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const double* p) {
			return _mm_loadu_pd (p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (double* p, Register a) {
			_mm_store_pd (p, a); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (double* p, Register a) {
			_mm_storeu_pd (p, a); 
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm_add_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register &a, const Register &b) {
			return _mm_add_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm_sub_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register &a, const Register &b) {
			return _mm_sub_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register &a, const Register &b) {
			return _mm_mul_pd(a,b);				
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return _mm_mul_sd(a,b);				
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &a, const Register &b) {
			return _mm_div_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm_div_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm_sqrt_pd(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register &a, const Register &b) {
			return _mm_sqrt_sd(a,b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm_min_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm_min_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm_max_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm_max_sd(a, b);
		}
		
	}; // SSETraits< double >

	template<>
	struct SSETraits<float> {
		
		typedef __m128 Register;         /**< @brief register type */
		static const unsigned int ne = 4; /**< @brief # of processed elements */
		static const unsigned int ns = 1; /**< @brief # of processed elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const float* p) {
			return _mm_load_ps (p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const float* p) {
			return _mm_loadu_ps (p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (float* p, Register a) {
			_mm_store_ps (p, a); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (float* p, Register a) {
			_mm_storeu_ps (p, a); 
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm_add_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register &a, const Register &b) {
			return _mm_add_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm_sub_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register &a, const Register &b) {
			return _mm_sub_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register &a, const Register &b) {
			return _mm_mul_ps(a, b);
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return _mm_mul_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &a, const Register &b) {
			return _mm_div_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm_div_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm_sqrt_ps(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register &a) {
			return _mm_sqrt_ss(a);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm_min_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm_min_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm_max_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm_max_ss(a, b);
		}
		
	}; // SSETraits<float>	

	template<>
	struct SSETraits< std::complex<float> > {
		
		typedef __m128 Register;         /**< @brief register type */
		static const unsigned int ne = 2; /**< @brief # of processed elements */
		static const unsigned int ns = 2; /**< @brief # of processed elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const std::complex<float>* p) {
			return _mm_load_ps ((float*)p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const std::complex<float>* p) {
			return _mm_loadu_ps ((float*)p); 
		}

		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (std::complex<float>* p, Register a) {
			_mm_store_ps ((float*)p, a); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (std::complex<float>* p, Register a) {
			_mm_storeu_ps ((float*)p, a); 
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register &a, const Register &b) {
			return _mm_add_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register &a, const Register &b) {
			return _mm_add_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register &a, const Register &b) {
			return _mm_sub_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register &a, const Register &b) {
			return _mm_sub_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register &ab, const Register &cd) {

			Register aa, bb, dc, x0, x1;  

			aa = _mm_moveldup_ps(ab);  
			bb = _mm_movehdup_ps(ab);  
			x0 = _mm_mul_ps(aa, cd);    //ac ad  
			dc = _mm_shuffle_ps(cd, cd, _MM_SHUFFLE(2,3,0,1));  
			x1 = _mm_mul_ps(bb, dc);    //bd bc  

			return _mm_addsub_ps(x0, x1);

		}
		
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return mulp(a, b);
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register &ab, const Register &cd) {

			Register aa, bb, dc, x0, x1, x2, x3;  

			bb = _mm_movehdup_ps (ab);   //  (   b1    b1    b0    b0) 
			aa = _mm_moveldup_ps (ab);   //  (   a1    a1    a0    a0) 

			dc = _mm_shuffle_ps(cd, cd, _MM_SHUFFLE(2,3,0,1));  
			                             //  (   c1    d1,   c0,   d0)  

			x0 = _mm_mul_ps (bb, cd);    //  ( b1d1  b1c1  b0d0  b0c0) 
			x1 = _mm_mul_ps (aa, dc);    //  ( a1c1  a1d1  a0c0  a0d0) 

			x2 = _mm_addsub_ps (x0, x1); //  ( b1d1  a1d1  b0d0  b0c0)
                                         //  (+a1c1 -a1d1 +a0c0 -a0d0)
			
			x1 = _mm_mul_ps (dc, dc);    //  ( c1c1  d1d1  c0c0  d0d0)
			x0 = _mm_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));
			                             //  ( d1d1  c1c1  d0d0  c0c0)

			x3 = _mm_add_ps (x1, x0);    //  ( c1c1  d1d1  c0c0  d0d0)
                                         //  (+d1d1  c1c1  d0d0  c0c0)

			x1 = _mm_div_ps (x2, x3);    //  ( b1d1+a1c1)
                                         //   /c1c1+d1d1)

			return _mm_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register &a, const Register &b) {
			return _mm_div_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register &a) {
			return _mm_sqrt_ps(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register &a) {
			return _mm_sqrt_ss(a);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register &a, const Register &b) {
			return _mm_min_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register &a, const Register &b) {
			return _mm_min_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register &a, const Register &b) {
			return _mm_max_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register &a, const Register &b) {
			return _mm_max_ss(a, b);
		}
		
	}; // SSETraits<float>	

#endif

} // namespace SSE 
