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

template<> struct SSETraits<float> {

    typedef float T;
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 8; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of processed elements */

    static inline Register loada (const T* p) {
        return _mm256_load_ps (p);
    }
    static inline Register loadu (const T* p) {
        return _mm256_loadu_ps (p);
    }
    static inline void stora (T* p, const Register& a) {
        _mm256_stream_ps (p, a);
    }
    static inline void storu (T* p, const Register& a) {
        _mm256_storeu_ps (p, a);
    }
    static inline Register addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }
    static inline Register adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }
    static inline Register subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }
    static inline Register subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }
    static inline Register mulp (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }
    static inline Register muls (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }
    static inline Register divp (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }
    static inline Register divs (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }
    static inline Register sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }
    static inline Register sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }
    static inline Register minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }
    static inline Register mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }
    static inline Register maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }
    static inline Register maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>

template<> struct SSETraits<double> {

    typedef double T;
    typedef __m256d Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */
    
    static inline Register loada (const T* p) {
        return _mm256_load_pd (p); 
    }
    static inline Register loadu (const T* p) {
        return _mm256_loadu_pd (p); 
    }
    static inline void stora (T* p, const Register& a) {
        _mm256_stream_pd (p, a);
    }
    static inline void storu (T* p, const Register& a) {
        _mm256_storeu_pd (p, a); 
    }
    static inline Register addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }
    static inline Register adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }
    static inline Register subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }
    static inline Register subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }
    static inline Register mulp (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }
    static inline Register muls (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }
    static inline Register divp (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }
    static inline Register divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }
    static inline Register sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }
    static inline Register sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }
    static inline Register minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }
    static inline Register mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }
    static inline Register maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }
    static inline Register maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< double >

template<>
struct SSETraits< std::complex<float> > {
    
    typedef std::complex<float> T;
    typedef float RT;
    
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of processed elements */
    

    static inline Register loada (const T* p) {
        return _mm256_load_ps (reinterpret_cast<const RT*>(p));
    }
    static inline Register loadu (const T* p) {
        return _mm256_loadu_ps (reinterpret_cast<const RT*>(p));
    }
    static inline void stora (T* p, const Register& a) {
        _mm256_stream_ps (reinterpret_cast<RT*>(p), a);
    }
    static inline void storu (T* p, const Register& a) {
        _mm256_storeu_ps (reinterpret_cast<RT*>(p), a); 
    }
    static inline Register addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }
    static inline Register adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }
    static inline Register subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }
    static inline Register subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }
    static inline Register mulp (const Register &a, const Register &b) {
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
    static inline Register muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }
    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    static inline Register change_sign (const Register & a) {
        if ((i0 | i1 | i2 | i3 | i4 | i5 | i6 | i7) == 0) return a;
        __m256 mask = constant8f<i0 ? 0x80000000 : 0, i1 ? 0x80000000 : 0, i2 ? 0x80000000 : 0, i3 ? 0x80000000 : 0,
                                 i4 ? 0x80000000 : 0, i5 ? 0x80000000 : 0, i6 ? 0x80000000 : 0, i7 ? 0x80000000 : 0> ();
        return _mm256_xor_ps(a, mask);
    }
    static inline Register divp (const Register &a, const Register &b) {
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
    static inline Register divs (const Register &a, const Register &b) {
        return divp(a, b);
    }
    static inline Register sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }
    static inline Register sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }
    static inline Register minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }
    static inline Register mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }
    static inline Register maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }
    static inline Register maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>    

template<>
struct SSETraits< std::complex<double> > {

    typedef std::complex<double> T;
    typedef double RT;
    typedef __m256d Register;
    static const unsigned int ne = 2;
    static const unsigned int ns = 2;

    static inline Register loada (const T* p) {
        return _mm256_load_pd (reinterpret_cast<const RT*>(p));
    }
    static inline Register loadu (const T* p) {
        return _mm256_loadu_pd (reinterpret_cast<const RT*>(p));
    }
    static inline void stora (T* p, const Register& a) {
        _mm256_stream_pd (reinterpret_cast<RT*>(p), a);
    }
    static inline void storu (T* p, const Register& a) {
        _mm256_storeu_pd (reinterpret_cast<RT*>(p), a);
    }
    static inline Register addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }
    static inline Register adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }
    static inline Register subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }
    static inline Register subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }
    static inline Register mulp (const Register &a, const Register &b) {
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
    static inline Register muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }
    template <int i0, int i1, int i2, int i3>
    static inline Register change_sign (const Register& a) {
        if ((i0 | i1 | i2 | i3) == 0)
            return a;
        __m256 mask = constant8f<0, i0 ? 0x80000000 : 0, 0, i1 ? 0x80000000 : 0, 0, i2 ? 0x80000000 : 0, 0, i3 ? 0x80000000 : 0> ();
        return _mm256_xor_pd(a, _mm256_castps_pd(mask));
    }
    static inline Register
    divp (const Register &a, const Register &b) {
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
    static inline Register divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }
    static inline Register sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }
    static inline Register sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }
    static inline Register minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }
    static inline Register mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }
    static inline Register maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }
    static inline Register maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< std::complex<double> >


#endif
