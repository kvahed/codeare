#ifndef __SSE2_TRAITS__
#define __SSE2_TRAITS__

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>


#include "SIMDTraits.hpp"

template <int i0, int i1, int i2, int i3>
static inline __m128i constant4i() {
    static const union {
        int     i[4];
        __m128i xmm;
    } u = {{i0,i1,i2,i3}};
    return u.xmm;
}

template<class T> struct SSETraits;


template<>
struct SSETraits< std::complex<double> > {
    
	typedef double RType;
	typedef std::complex<double> Type;
    typedef __m128d Register;         /**< @brief register type */
    static const unsigned int ne = 1; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of sub elements */
	
    /**
     * @brief     SSE2 load packed aligned
     */
    static inline Register 
    loada (const Type* p) {
        return _mm_load_pd ((const RType*)p);
    }
    
    /**
     * @brief     SSE2 load packed unaligned
     */
    static inline Register
    loadu (const Type* p) {
        return _mm_loadu_pd ((const RType*)p);
    }
    
    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
        _mm_stream_pd ((RType*)p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
        _mm_storeu_pd ((RType*)p, a); 
    }

    /**
     * @brief     SSE2 packed addition
     */
    static inline Register
    addp (const Register& a, const Register& b) {
        return _mm_add_pd(a, b);
    }
	
    /**
     * @brief     SSE2 single addition
     */
    static inline Register 
    adds (const Register& a, const Register& b) {
        return _mm_add_sd(a, b);
    }
	
    /**
     * @brief     SSE2 packed subtraction
     */
    static inline Register 
    subp (const Register& a, const Register& b) {
        return _mm_sub_pd(a, b);
    }
	
    /**
     * @brief     SSE2 single subtraction
     */
    static inline Register 
    subs (const Register& a, const Register& b) {
        return _mm_sub_sd(a, b);
    }
	
    template <int i0, int i1>
    static inline Register change_sign(const Register& a) {
        if ((i0 | i1) == 0) return a;
        __m128i mask = constant4i<0, i0 ? 0x80000000 : 0, 0, i1 ? 0x80000000 : 0> ();
        return  _mm_xor_pd(a, _mm_castsi128_pd(mask));     // flip sign bits
    }

    /**
     * @brief     SSE3 packed multiplication
     */
    static inline Register 
    mulp (const Register& a, const Register& b) {
        
        Register b_flip = _mm_shuffle_pd(b,b,1);      // Swap b.re and b.im
        Register a_im   = _mm_shuffle_pd(a,a,3);      // Imag. part of a in both
        Register a_re   = _mm_shuffle_pd(a,a,0);      // Real part of a in both
        Register aib    = _mm_mul_pd(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
#ifdef __FMA__      // FMA3
        return  _mm_fmaddsub_pd(a_re, b, aib);       // a_re * b +/- aib
#elif defined (__FMA4__)  // FMA4
        return  _mm_maddsub_pd (a_re, b, aib);       // a_re * b +/- aib
#elif  INSTRSET >= 3  // SSE3
        Register arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        return _mm_addsub_pd(arb, aib);              // subtract/add
#else
        Register arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        Register aib_m  = change_sign<1,0>(Register(aib)); // change sign of low part
        return _mm_add_pd(arb, aib_m);               // add
#endif
        
    }
	
    /**
     * @brief     SSE2 single multiplication
     */
    static inline Register 
    muls (const Register& a, const Register& b) {
        return mulp(a, b);
    }
	/*
    template <int i0, int i1>
    static inline Register change_sign(Register const & a) {
        if ((i0 | i1) == 0) return a;
        __m128i mask = constant4i<0, i0 ? 0x80000000 : 0, 0, i1 ? 0x80000000 : 0> ();
        return  _mm_xor_pd(a, _mm_castsi128_pd(mask));     // flip sign bits
    }
    */
    /**
     * @brief     SSE2 packed division
     */
    static inline Register 
    divp (const Register& a, const Register& b) {
        
        Register a_re   = _mm_shuffle_pd(a,a,0);      // Real part of a in both
        Register arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        Register b_flip = _mm_shuffle_pd(b,b,1);      // Swap b.re and b.im
        Register a_im   = _mm_shuffle_pd(a,a,3);      // Imag. part of a in both
#ifdef __FMA__      // FMA3
        Register n      = _mm_fmsubadd_pd(a_im, b_flip, arb);
#elif defined (__FMA4__)  // FMA4
        Register n      = _mm_msubadd_pd (a_im, b_flip, arb);
#else
        Register aib    = _mm_mul_pd(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
        Register arbm   = change_sign<0,1>(Register(arb));
        Register n      = _mm_add_pd(arbm, aib);      // arbm + aib
#endif  // FMA
        Register bb     = _mm_mul_pd(b, b);           // (b.re*b.re, b.im*b.im)
#if INSTRSET >= 3  // SSE3
        Register bsq    = _mm_hadd_pd(bb,bb);         // (b.re*b.re + b.im*b.im) dublicated
#else
        Register bb1    = _mm_shuffle_pd(bb,bb,1);
        Register bsq    = _mm_add_pd(bb,bb1);
#endif // SSE3
        return           _mm_div_pd(n, bsq);         // n / bsq
        
    }
	
    /**
     * @brief     SSE2 single division
     */
    static inline Register 
    divs (const Register& a, const Register& b) {
        return _mm_div_sd(a, b);
    }
    
    /**
     * @brief     SSE2 packed SQRT
     */
    static inline Register 
    sqrtp (const Register& a) {
        return _mm_sqrt_pd(a);
    }
	
    /**
     * @brief     SSE2 single SQRT
     */
    static inline Register 
    sqrts (const Register& a, const Register& b) {
        return _mm_sqrt_sd(a,b);
    }
	
    /**
     * @brief     SSE2 packed comparison
     */
    static inline Register 
    minp (const Register& a, const Register& b) {
        return _mm_min_pd(a, b);
    }
	
    /**
     * @brief     SSE2 single comparison
     */
    static inline Register 
    mins (const Register& a, const Register& b) {
        return _mm_min_sd(a, b);
    }
	
    /**
     * @brief     SSE2 packed comparison
     */
    static inline Register 
    maxp (const Register& a, const Register& b) {
        return _mm_max_pd(a, b);
    }
		
    /**
     * @brief     SSE2 single comparison
     */
    static inline Register 
    maxs (const Register& a, const Register& b) {
			return _mm_max_sd(a, b);
    }
	
}; // SSETraits< std::complex<double> >

template<>
	struct SSETraits< double > {
		
		typedef double Type;
		typedef __m128d Register;         /**< @brief register type */
		static const unsigned int ne = 2; /**< @brief # of processed elements */
		static const unsigned int ns = 1; /**< @brief # of sub elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const Type* p) {
			return _mm_load_pd (p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const Type* p) {
			return _mm_loadu_pd (p); 
		}

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
        _mm_stream_pd (p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
        _mm_storeu_pd (p, a); 
    }

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register& a, const Register& b) {
			return _mm_add_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register& a, const Register& b) {
			return _mm_add_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register& a, const Register& b) {
			return _mm_sub_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register& a, const Register& b) {
			return _mm_sub_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register& a, const Register& b) {
			return _mm_mul_pd(a,b);				
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register& a, const Register& b) {
			return _mm_mul_sd(a,b);				
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register& a, const Register& b) {
			return _mm_div_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register& a, const Register& b) {
			return _mm_div_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register& a) {
			return _mm_sqrt_pd(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register& a, const Register& b) {
			return _mm_sqrt_sd(a,b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register& a, const Register& b) {
			return _mm_min_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register& a, const Register& b) {
			return _mm_min_sd(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register& a, const Register& b) {
			return _mm_max_pd(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register& a, const Register& b) {
			return _mm_max_sd(a, b);
		}
		
	}; // SSETraits< double >

	template<>
	struct SSETraits<float> {
		
		typedef float Type;
		typedef __m128 Register;         /**< @brief register type */
		static const unsigned int ne = 4; /**< @brief # of processed elements */
		static const unsigned int ns = 1; /**< @brief # of processed elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const Type* p) {
			return _mm_load_ps (p); 
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const Type* p) {
			return _mm_loadu_ps (p); 
		}

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
        _mm_stream_ps (p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
        _mm_storeu_ps (p, a); 
    }


		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register& a, const Register& b) {
			return _mm_add_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register& a, const Register& b) {
			return _mm_add_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register& a, const Register& b) {
			return _mm_sub_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register& a, const Register& b) {
			return _mm_sub_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register& a, const Register& b) {
			return _mm_mul_ps(a, b);
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register& a, const Register& b) {
			return _mm_mul_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register& a, const Register& b) {
			return _mm_div_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register& a, const Register& b) {
			return _mm_div_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register& a) {
			return _mm_sqrt_ps(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register& a) {
			return _mm_sqrt_ss(a);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register& a, const Register& b) {
			return _mm_min_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register& a, const Register& b) {
			return _mm_min_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register& a, const Register& b) {
			return _mm_max_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register& a, const Register& b) {
			return _mm_max_ss(a, b);
		}
		
	}; // SSETraits<float>	

	template<>
	struct SSETraits< std::complex<float> > {
		
		typedef std::complex<float> Type;
		typedef float RType;
		typedef __m128 Register;         /**< @brief register type */
		static const unsigned int ne = 2; /**< @brief # of processed elements */
		static const unsigned int ns = 2; /**< @brief # of processed elements */
		
		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline Register 
		loada (const Type* p) {
			return _mm_load_ps ((const RType*)p);
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline Register
		loadu (const Type* p) {
			return _mm_loadu_ps ((const RType*)p);
		}


		/**
		 * @brief     SSE2 load packed aligned
		 */
		static inline void
		stora (Type* p, Register a) {
			_mm_stream_ps ((RType*)p, a);
		}

		/**
		 * @brief     SSE2 load packed unaligned
		 */
		static inline void
		storu (Type* p, Register a) {
			_mm_storeu_ps ((RType*)p, a);
		}

		/**
		 * @brief     SSE2 packed addition
		 */
		static inline Register
		addp (const Register& a, const Register& b) {
			return _mm_add_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single addition
		 */
		static inline Register 
		adds (const Register& a, const Register& b) {
			return _mm_add_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed subtraction
		 */
		static inline Register 
		subp (const Register& a, const Register& b) {
			return _mm_sub_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single subtraction
		 */
		static inline Register 
		subs (const Register& a, const Register& b) {
			return _mm_sub_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed multiplication
		 */
		static inline Register 
		mulp (const Register& ab, const Register& cd) {

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
		muls (const Register& a, const Register& b) {
			return mulp(a, b);
		}
		
		/**
		 * @brief     SSE2 packed division
		 */
		static inline Register 
		divp (const Register& ab, const Register& cd) {

			Register aa, bb, dc, x0, x1, x2, x3;  

			bb = _mm_movehdup_ps (ab);
			aa = _mm_moveldup_ps (ab);
			dc = _mm_shuffle_ps(cd, cd, _MM_SHUFFLE(2,3,0,1));  
			x0 = _mm_mul_ps (bb, cd);
			x1 = _mm_mul_ps (aa, dc);
			x2 = _mm_addsub_ps (x0, x1);
			x1 = _mm_mul_ps (dc, dc);
			x0 = _mm_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));
			x3 = _mm_add_ps (x1, x0);
			x1 = _mm_div_ps (x2, x3);

			return _mm_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));

		}
		
		/**
		 * @brief     SSE2 single division
		 */
		static inline Register 
		divs (const Register& a, const Register& b) {
			return _mm_div_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed SQRT
		 */
		static inline Register 
		sqrtp (const Register& a) {
			return _mm_sqrt_ps(a);
		}
		
		/**
		 * @brief     SSE2 single SQRT
		 */
		static inline Register 
		sqrts (const Register& a) {
			return _mm_sqrt_ss(a);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		minp (const Register& a, const Register& b) {
			return _mm_min_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		mins (const Register& a, const Register& b) {
			return _mm_min_ss(a, b);
		}
		
		/**
		 * @brief     SSE2 packed comparison
		 */
		static inline Register 
		maxp (const Register& a, const Register& b) {
			return _mm_max_ps(a, b);
		}
		
		/**
		 * @brief     SSE2 single comparison
		 */
		static inline Register 
		maxs (const Register& a, const Register& b) {
			return _mm_max_ss(a, b);
		}
		
	}; // SSETraits<float>	

#endif
