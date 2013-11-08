#ifndef __AVX_TRAITS__
#define __AVX_TRAITS__

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>

template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
static inline __m256 constant8f() {
    static const union {
        int     i[8];
        __m256  ymm;
    } u = {{i0,i1,i2,i3,i4,i5,i6,i7}};
    return u.ymm;
}

template<class T> struct SSETraits;

template<>
struct SSETraits< std::complex<double> > {

    typedef std::complex<double> type;
    typedef __m256d Register;         /**< @brief register type */
    static const unsigned int ne = 2; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_pd ((double*)p);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_pd ((double*)p);
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loadoa (const type* p) {
        return _mm256_load_pd ((double*)p);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadou (const type* p) {
        return _mm256_loadu_pd ((double*)p);
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_stream_pd ((double*)p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_pd ((double*)p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     SSE3 packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {
        
        __m256d b_flip = _mm256_shuffle_pd(b,b,5); // Swap b.re and b.im
        __m256d a_im = _mm256_shuffle_pd(a,a,0xF); // Imag part of a in both
        __m256d a_re = _mm256_shuffle_pd(a,a,0); // Real part of a in both
        __m256d aib = _mm256_mul_pd(a_im, b_flip); // (a.im*b.im, a.im*b.re)
#ifdef __FMA__ // FMA3
        return _mm256_fmaddsub_pd(a_re, b, aib); // a_re * b +/- aib
#elif defined (__FMA4__) // FMA4
        return _mm256_maddsub_pd (a_re, b, aib); // a_re * b +/- aib
#else
        __m256d arb = _mm256_mul_pd(a_re, b); // (a.re*b.re, a.re*b.im)
        return _mm256_addsub_pd(arb, aib); // subtract/add
#endif // FMA
        
    }
    
    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }

    template <int i0, int i1, int i2, int i3>
    static inline Register
    change_sign (const Register& a) {
        if ((i0 | i1 | i2 | i3) == 0)
            return a;
        __m256 mask = constant8f<0, i0 ? 0x80000000 : 0, 0, i1 ? 0x80000000 : 0, 0, i2 ? 0x80000000 : 0, 0, i3 ? 0x80000000 : 0> ();
        return _mm256_xor_pd(a, _mm256_castps_pd(mask));
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        
        // The following code is made similar to the operator * to enable common
        // subexpression elimination in code that contains both operator * and
        // operator / where one or both operands are the same
        __m256d a_re = _mm256_shuffle_pd(a,a,0); // Real part of a in both
        __m256d arb = _mm256_mul_pd(a_re, b); // (a.re*b.re, a.re*b.im)
        __m256d b_flip = _mm256_shuffle_pd(b,b,5); // Swap b.re and b.im
        __m256d a_im = _mm256_shuffle_pd(a,a,0xF); // Imag part of a in both
#ifdef __FMA__ // FMA3
        __m256d n = _mm256_fmsubadd_pd(a_im, b_flip, arb);
#elif defined (__FMA4__) // FMA4
        __m256d n = _mm256_msubadd_pd (a_im, b_flip, arb);
#else
        __m256d aib = _mm256_mul_pd(a_im, b_flip); // (a.im*b.im, a.im*b.re)
        __m256d arbm = change_sign<0,1,0,1>(Register(arb));
        __m256d n = _mm256_add_pd(arbm, aib); // arbm + aib
#endif // FMA
        __m256d bb = _mm256_mul_pd(b, b); // (b.re*b.re, b.im*b.im)
        __m256d bsq = _mm256_hadd_pd(bb,bb); // (b.re*b.re + b.im*b.im) dublicated
        return _mm256_div_pd(n, bsq); // n / 
        
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< std::complex<double> >

template<>
struct SSETraits<double> {

    typedef __m256d Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const double* p) {
        return _mm256_load_pd (p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const double* p) {
        return _mm256_loadu_pd (p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (double* p, Register a) {
        _mm256_stream_pd (p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (double* p, Register a) {
        _mm256_storeu_pd (p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< double >

template<>
struct SSETraits<float> {

    typedef float type;
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 8; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of processed elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_ps (p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_ps (p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_stream_ps (p, a);
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_ps (p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>    

template<>
struct SSETraits< std::complex<float> > {
    
    typedef std::complex<float> type;
    typedef float rtype;
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of processed elements */
    
    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_ps ((rtype*)p);
    }
    
    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_ps ((rtype*)p);
    }
    
    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_stream_ps ((rtype*)p, a);
    }
    
    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_ps ((rtype*)p, a);
    }
    
    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {

    __m256 b_flip = _mm256_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
    __m256 a_im   = _mm256_shuffle_ps(a,a,0xF5);   // Imag part of a in both
    __m256 a_re   = _mm256_shuffle_ps(a,a,0xA0);   // Real part of a in both
    __m256 aib    = _mm256_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
#ifdef __FMA__      // FMA3
    return  _mm256_fmaddsub_ps(a_re, b, aib);      // a_re * b +/- aib
#elif defined (__FMA4__)  // FMA4
    return  _mm256_maddsub_ps (a_re, b, aib);      // a_re * b +/- aib
#else
    __m256 arb    = _mm256_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
    return          _mm256_addsub_ps(arb, aib);    // subtract/add
#endif  // FMA

    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }


    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    static inline Register change_sign (const Register & a) {
        if ((i0 | i1 | i2 | i3 | i4 | i5 | i6 | i7) == 0) return a;
        __m256 mask = constant8f<i0 ? 0x80000000 : 0, i1 ? 0x80000000 : 0, i2 ? 0x80000000 : 0, i3 ? 0x80000000 : 0,
                                 i4 ? 0x80000000 : 0, i5 ? 0x80000000 : 0, i6 ? 0x80000000 : 0, i7 ? 0x80000000 : 0> ();
        return _mm256_xor_ps(a, mask);
    }


    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        
        __m256 a_re   = _mm256_shuffle_ps(a,a,0xA0);   // Real part of a in both
        __m256 arb    = _mm256_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
        __m256 b_flip = _mm256_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
        __m256 a_im   = _mm256_shuffle_ps(a,a,0xF5);   // Imag part of a in both
#ifdef __FMA__      // FMA3
        __m256 n      = _mm256_fmsubadd_ps(a_im, b_flip, arb); 
#elif defined (__FMA4__)  // FMA4
        __m256 n      = _mm256_msubadd_ps (a_im, b_flip, arb);
#else
        __m256 aib    = _mm256_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
        __m256 arbm   = change_sign<0,1,0,1,0,1,0,1>(Register(arb)); // (a.re*b.re,-a.re*b.im)
        __m256 n      = _mm256_add_ps(arbm, aib);      // arbm + aib
#endif  // FMA
        __m256 bb     = _mm256_mul_ps(b, b);           // (b.re*b.re, b.im*b.im)
        __m256 bb2    = _mm256_shuffle_ps(bb,bb,0xB1); // Swap bb.re and bb.im
        __m256 bsq    = _mm256_add_ps(bb,bb2);         // (b.re*b.re + b.im*b.im) dublicated
        return          _mm256_div_ps(n, bsq);         // n / bsq

    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return divp(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>    


#endif
