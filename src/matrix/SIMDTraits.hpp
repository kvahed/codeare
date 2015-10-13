/*
 * SIMDTraits.hpp
 *
 *  Created on: Sep 2, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_SIMDTRAITS_HPP_
#define SRC_MATRIX_SIMDTRAITS_HPP_

#include "Vector.hpp"
#include "TypeTraits.hpp"
#include <algorithm>

template<class T> struct VecTraits;

#if defined (__AVX__)
#include <immintrin.h>

template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
static inline __m256 constant8f() {
    static const union {
        int     i[8];
        __m256  ymm;
    } u = {{i0,i1,i2,i3,i4,i5,i6,i7}};
    return u.ymm;
}

template<> struct VecTraits<float> {
    typedef __m256 reg_type;
    static const int stride = 8;
    inline static reg_type setp (const float& s) { return _mm256_set1_ps(s); }
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_ps(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm256_mul_ps(a, b);}
    inline static reg_type divides (const reg_type& a, const reg_type& b) {return _mm256_div_ps(a, b);}
    inline static reg_type conjugate (const reg_type& a) { return a; }
};

template<> struct VecTraits<double> {
    typedef __m256d reg_type;
    static const int stride = 4;
    inline static reg_type setp (const double& s) { return _mm256_set1_pd(s); }
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_pd(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm256_mul_pd(a, b);}
    inline static reg_type divides (const reg_type& a, const reg_type& b) {return _mm256_div_pd(a, b);}
    inline static reg_type conjugate (const reg_type& a) { return a; }
};

template<> struct VecTraits<cxfl> {
    typedef __m256 reg_type;
    static const int stride = 4;
    inline static reg_type setp (const cxfl& s) {
        return _mm256_set_ps(real(s),imag(s),real(s),imag(s),real(s),imag(s),real(s),imag(s)); }
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_ps(a, b);}
    inline static reg_type multiplies (const reg_type &a, const reg_type &b) {
    	reg_type b_flip = _mm256_shuffle_ps(b,b,0xB1); // Swap b.re and b.im
    	reg_type a_im   = _mm256_shuffle_ps(a,a,0xF5); // Imag part of a in both
    	reg_type a_re   = _mm256_shuffle_ps(a,a,0xA0); // Real part of a in both
    	reg_type aib    = _mm256_mul_ps(a_im, b_flip); // (a.im*b.im, a.im*b.re)
#ifdef __FMA__           // FMA3
    	return  _mm256_fmaddsub_ps(a_re, b, aib);      // a_re * b +/- aib
#elif defined (__FMA4__) // FMA4
    	return  _mm256_maddsub_ps (a_re, b, aib);      // a_re * b +/- aib
#else
    	reg_type arb    = _mm256_mul_ps(a_re, b);      // (a.re*b.re, a.re*b.im)
    	return          _mm256_addsub_ps(arb, aib);    // subtract/add
#endif                   // FMA
    }
    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    inline static reg_type change_sign (const reg_type & a) {
        if ((i0 | i1 | i2 | i3 | i4 | i5 | i6 | i7) == 0) return a;
        reg_type mask = constant8f<i0 ? 0x80000000 : 0, i1 ? 0x80000000 : 0, i2 ? 0x80000000 : 0, i3 ? 0x80000000 : 0,
                                 i4 ? 0x80000000 : 0, i5 ? 0x80000000 : 0, i6 ? 0x80000000 : 0, i7 ? 0x80000000 : 0> ();
        return _mm256_xor_ps(a, mask);
    }
    inline static reg_type divides (const reg_type &a, const reg_type &b) {
        reg_type a_re   = _mm256_shuffle_ps(a,a,0xA0);   // Real part of a in both
        reg_type arb    = _mm256_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
        reg_type b_flip = _mm256_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
        reg_type a_im   = _mm256_shuffle_ps(a,a,0xF5);   // Imag part of a in both
#ifdef __FMA__      // FMA3
        reg_type n      = _mm256_fmsubadd_ps(a_im, b_flip, arb); 
#elif defined (__FMA4__)  // FMA4
        reg_type n      = _mm256_msubadd_ps (a_im, b_flip, arb);
#else
        reg_type aib    = _mm256_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
        reg_type arbm   = change_sign<0,1,0,1,0,1,0,1>(arb); // (a.re*b.re,-a.re*b.im)
        reg_type n      = _mm256_add_ps(arbm, aib);      // arbm + aib
#endif  // FMA
        reg_type bb     = _mm256_mul_ps(b, b);           // (b.re*b.re, b.im*b.im)
        reg_type bb2    = _mm256_shuffle_ps(bb,bb,0xB1); // Swap bb.re and bb.im
        reg_type bsq    = _mm256_add_ps(bb,bb2);         // (b.re*b.re + b.im*b.im) dublicated
        return          _mm256_div_ps(n, bsq);         // n / bsq
    }
  inline static reg_type conjugate (const reg_type& a) {
    return change_sign<0,1,0,1,0,1,0,1>(a);
  }
};

template<> struct VecTraits<cxdb> {
    typedef __m256d reg_type;
    static const int stride = 2;
    inline static reg_type setp (const cxdb& s) {
        return _mm256_set_pd(real(s),imag(s),real(s),imag(s)); }
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_pd(a, b);}
    inline static reg_type multiplies (const reg_type &a, const reg_type &b) {
        reg_type b_flip = _mm256_shuffle_pd(b,b,5);    // Swap b.re and b.im
        reg_type a_im   = _mm256_shuffle_pd(a,a,0xF);  // Imag part of a in both
        reg_type a_re   = _mm256_shuffle_pd(a,a,0);    // Real part of a in both
        reg_type aib    = _mm256_mul_pd(a_im, b_flip); // (a.im*b.im, a.im*b.re)
#ifdef __FMA__           // FMA3
        return _mm256_fmaddsub_pd(a_re, b, aib);       // a_re * b +/- aib
#elif defined (__FMA4__) // FMA4
        return _mm256_maddsub_pd (a_re, b, aib);       // a_re * b +/- aib
#else
        reg_type arb    = _mm256_mul_pd(a_re, b);      // (a.re*b.re, a.re*b.im)
        return _mm256_addsub_pd(arb, aib);             // subtract/add
#endif                   // FMA
    }
    template <int i0, int i1, int i2, int i3>
    inline static reg_type change_sign (const reg_type& a) {
        if ((i0 | i1 | i2 | i3) == 0)
            return a;
        __m256 mask = constant8f<0, i0 ? 0x80000000 : 0, 0, i1 ? 0x80000000 : 0, 0, i2 ? 0x80000000 : 0, 0, i3 ? 0x80000000 : 0> ();
        return _mm256_xor_pd(a, _mm256_castps_pd(mask));
    }
    inline static reg_type divides (const reg_type &a, const reg_type &b) {
        reg_type a_re   = _mm256_shuffle_pd(a,a,0); // Real part of a in both
        reg_type arb    = _mm256_mul_pd(a_re, b); // (a.re*b.re, a.re*b.im)
        reg_type b_flip = _mm256_shuffle_pd(b,b,5); // Swap b.re and b.im
        reg_type a_im   = _mm256_shuffle_pd(a,a,0xF); // Imag part of a in both
#ifdef __FMA__ // FMA3
        reg_type n      = _mm256_fmsubadd_pd(a_im, b_flip, arb);
#elif defined (__FMA4__) // FMA4
        reg_type n      = _mm256_msubadd_pd (a_im, b_flip, arb);
#else
        reg_type aib    = _mm256_mul_pd(a_im, b_flip); // (a.im*b.im, a.im*b.re)
        reg_type arbm   = change_sign<0,1,0,1>(reg_type(arb));
        reg_type n      = _mm256_add_pd(arbm, aib); // arbm + aib
#endif // FMA
        reg_type bb     = _mm256_mul_pd(b, b); // (b.re*b.re, b.im*b.im)
        reg_type bsq    = _mm256_hadd_pd(bb,bb); // (b.re*b.re + b.im*b.im) dublicated
        return _mm256_div_pd(n, bsq); // n /
    }
  inline static reg_type conjugate (const reg_type& a) {
    return change_sign<0,1,0,1>(a);
  }
};
#elif defined (__SSE2__)
#include <emmintrin.h>

template <int i0, int i1, int i2, int i3>
inline static __m128i constant4i() {
    static const union {
        int     i[4];
        __m128i xmm;
    } u = {{i0,i1,i2,i3}};
    return u.xmm;
}

template <int i0, int i1, int i2, int i3>
inline static __m128 change_sign(__m128 const & a) {
    if ((i0 | i1 | i2 | i3) == 0) return a;
    __m128i mask = constant4i<i0 ? 0x80000000 : 0, i1 ? 0x80000000 : 0, i2 ? 0x80000000 : 0, i3 ? 0x80000000 : 0>();
    return  _mm_xor_ps(a, _mm_castsi128_ps(mask));     // flip sign bits
}

template<> struct VecTraits<float> {
    typedef __m128 reg_type;
    static const int stride = 4;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm_sub_ps(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm_mul_ps(a, b);}
    inline static reg_type divides (const reg_type& a, const reg_type& b) {return _mm_div_ps(a, b);}
    inline static reg_type conjugate (const reg_type& a) { return a; }

};

template<> struct VecTraits<double> {
    typedef __m128d reg_type;
    static const int stride = 2;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm_sub_pd(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm_mul_pd(a, b);}
    inline static reg_type divides (const reg_type& a, const reg_type& b) {return _mm_div_pd(a, b);}
    inline static reg_type conjugate (const reg_type& a) { return a; }
};
template<> struct VecTraits<cxfl> {
    typedef __m128 reg_type;
    static const int stride = 2;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm_sub_ps(a, b);}
    template <int i0, int i1, int i2, int i3>
    inline static reg_type change_sign(reg_type const & a) {
        if ((i0 | i1 | i2 | i3) == 0) return a;
        __m128i mask = constant4i<i0 ? 0x80000000 : 0, i1 ? 0x80000000 : 0, i2 ? 0x80000000 : 0, i3 ? 0x80000000 : 0>();
        return  _mm_xor_ps(a, _mm_castsi128_ps(mask));     // flip sign bits
    }
    inline static reg_type multiplies (reg_type const & a, reg_type const & b) {
        __m128 b_flip = _mm_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
        __m128 a_im   = _mm_shuffle_ps(a,a,0xF5);   // Imag. part of a in both
        __m128 a_re   = _mm_shuffle_ps(a,a,0xA0);   // Real part of a in both
        __m128 aib    = _mm_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
    #if defined (__FMA__)
        return  _mm_fmaddsub_ps(a_re, b, aib);      // a_re * b +/- aib
    #elif defined (__FMA4__)
        return  _mm_maddsub_ps (a_re, b, aib);      // a_re * b +/- aib
    #elif defined (__SSSE3__)
        __m128 arb    = _mm_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
        return _mm_addsub_ps(arb, aib);             // subtract/add
    #else
        __m128 arb     = _mm_mul_ps(a_re, b);       // (a.re*b.re, a.re*b.im)
        __m128 aib_m   = change_sign<1,0,1,0>(aib); // change sign of low part
        return _mm_add_ps(arb, aib_m);              // add
    #endif
    }
    inline static reg_type divides (reg_type const & a, reg_type const & b) {
        __m128 a_re   = _mm_shuffle_ps(a,a,0xA0);   // Real part of a in both
        __m128 arb    = _mm_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
        __m128 b_flip = _mm_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
        __m128 a_im   = _mm_shuffle_ps(a,a,0xF5);   // Imag. part of a in both
    #if defined (__FMA__)
        __m128 n      = _mm_fmsubadd_ps(a_im, b_flip, arb);
    #elif defined (__FMA4__)
        __m128 n      = _mm_msubadd_ps (a_im, b_flip, arb);
    #else
        __m128 aib    = _mm_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
        __m128 arbm   = change_sign<0,1,0,1>(arb); // change sign of high part
        __m128 n      = _mm_add_ps(arbm, aib);      // arbm + aib
    #endif  // FMA
        __m128 bb     = _mm_mul_ps(b, b);           // (b.re*b.re, b.im*b.im)
        __m128 bb1    = _mm_shuffle_ps(bb,bb,0xB1); // Swap bb.re and bb.im
        __m128 bb2    = _mm_add_ps(bb,bb1);         // add pairwise into both positions
        return          _mm_div_ps(n, bb2);         // n / bb3
    }
  inline static reg_type conjugate (const reg_type& a) {
    return change_sign<0,1,0,1>(a);
  }
};
template<> struct VecTraits<cxdb> {
    typedef __m128d reg_type;
    static const int stride = 2;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm_sub_pd(a, b);}
    inline static reg_type multiplies (reg_type const & a, reg_type const & b) {
        reg_type b_flip = _mm_shuffle_pd(b,b,1);      // Swap b.re and b.im
        reg_type a_im   = _mm_shuffle_pd(a,a,3);      // Imag. part of a in both
        reg_type a_re   = _mm_shuffle_pd(a,a,0);      // Real part of a in both
        reg_type aib    = _mm_mul_pd(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
    #if defined (__FMA__)      // FMA3
        return  _mm_fmaddsub_pd(a_re, b, aib);       // a_re * b +/- aib
    #elif defined (__FMA4__)  // FMA4
        return  _mm_maddsub_pd (a_re, b, aib);       // a_re * b +/- aib
    #elif  defined (__SSSE3__)  // SSE3
        reg_type arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        return _mm_addsub_pd(arb, aib);              // subtract/add
    #else
        reg_type arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        reg_type aib_m  = change_sign<1,0>(aib); // change sign of low part
        return _mm_add_pd(arb, aib_m);               // add
    #endif
    }
    inline static reg_type divides (reg_type const & a, reg_type const & b) {
        reg_type a_re   = _mm_shuffle_pd(a,a,0);      // Real part of a in both
        reg_type arb    = _mm_mul_pd(a_re, b);        // (a.re*b.re, a.re*b.im)
        reg_type b_flip = _mm_shuffle_pd(b,b,1);      // Swap b.re and b.im
        reg_type a_im   = _mm_shuffle_pd(a,a,3);      // Imag. part of a in both
#ifdef __FMA__      // FMA3
        reg_type n      = _mm_fmsubadd_pd(a_im, b_flip, arb);
#elif defined (__FMA4__)  // FMA4
        reg_type n      = _mm_msubadd_pd (a_im, b_flip, arb);
        #else
        reg_type aib    = _mm_mul_pd(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
        reg_type arbm   = change_sign<0,1>(Vec2d(arb));
        reg_type n      = _mm_add_pd(arbm, aib);      // arbm + aib
#endif  // FMA
        reg_type bb     = _mm_mul_pd(b, b);           // (b.re*b.re, b.im*b.im)
#if INSTRSET >= 3  // SSE3
        reg_type bsq    = _mm_hadd_pd(bb,bb);         // (b.re*b.re + b.im*b.im) dublicated
        #else
        reg_type bb1    = _mm_shuffle_pd(bb,bb,1);
        reg_type bsq    = _mm_add_pd(bb,bb1);
#endif // SSE3
        return           _mm_div_pd(n, bsq);         // n / bsq
    }
  inline static reg_type conjugate (const reg_type& a) {
    return change_sign<0,1>(a);
  }
};
#endif

namespace codeare {
	template<class T> class plus {
	public:
		typedef typename VecTraits<T>::reg_type reg_type;
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        inline static reg_type packed (const reg_type& a, const reg_type& b) {
			return VecTraits<T>::plus(a, b);
		}
		inline T operator() (const T& x, const T& y) const {
			return std::plus<T>()(x, y);
		}
	};
	template<class T> class minus {
	public:
		typedef typename VecTraits<T>::reg_type reg_type;
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
		inline static reg_type packed (const reg_type& a, const reg_type& b) {
			return VecTraits<T>::minus(a, b);
		}
		inline T operator() (const T& x, const T& y) const {
			return std::minus<T>()(x, y);
		}
	};
	template<class T> class multiplies {
	public:
		typedef typename VecTraits<T>::reg_type reg_type;
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        inline static reg_type packed (const reg_type& a, const reg_type& b) {
			return VecTraits<T>::multiplies(a, b);
		}
		inline T operator() (const T& x, const T& y) const {
			return std::multiplies<T>()(x, y);
		}
	};
	template<class T> class divides {
	public:
		typedef typename VecTraits<T>::reg_type reg_type;
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
		inline static reg_type packed (const reg_type& a, const reg_type& b) {
			return VecTraits<T>::divides(a, b);
		}
		inline T operator() (const T& x, const T& y) const {
			return std::divides<T>()(x, y);
		}
	};
	template<class T> class conjugate {
	public:
		typedef typename VecTraits<T>::reg_type reg_type;
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
		inline static reg_type packed (const reg_type& a) {
			return VecTraits<T>::conjugate(a);
		}
		T operator() (const T& x) const;
	};
	template<> inline float conjugate<float>::operator() (const float& f) const { return f; }
	template<> inline double conjugate<double>::operator() (const double& d) const { return d; }
	template<> inline cxfl conjugate<cxfl>::operator() (const cxfl& cf) const { return std::conj(cf); }
	template<> inline cxdb conjugate<cxdb>::operator() (const cxdb& cd) const { return std::conj(cd); }
}


template<class T, class Op> 
inline static void Vec (const Vector<T>& a, const Vector<T>& b, Vector<T>& c, const Op& op) {
    typedef typename VecTraits<T>::reg_type reg_type;
    const reg_type* va = (const reg_type*) &a[0];
    const reg_type* vb = (const reg_type*) &b[0];
    reg_type* vc = (reg_type*) &c[0];
    int simd_n = std::floor(a.size()/VecTraits<T>::stride);
    size_t start_r = simd_n*VecTraits<T>::stride;
    for (size_t i = 0; i < simd_n; ++i)
        vc[i] = op.packed(va[i], vb[i]);
    std::transform(a.begin()+start_r, a.end(), b.begin()+start_r, c.begin()+start_r, op);
}

template<class T, class Op> 
inline static void Vec (const Vector<T>& a, Vector<T>& c, const Op& op) {
    typedef typename VecTraits<T>::reg_type reg_type;
    const reg_type* va = (const reg_type*) &a[0];
    reg_type* vc = (reg_type*) &c[0];
    size_t simd_n = std::floor(a.size()/VecTraits<T>::stride),
      start_r = simd_n*VecTraits<T>::stride;
    for (size_t i = 0; i < simd_n; ++i)
        vc[i] = op.packed(va[i]);
    std::transform(a.begin()+start_r, a.end(), c.begin()+start_r, op);
}

template<class T, class Op> 
inline static void Vec (const Vector<T>& a, const T& b, Vector<T>& c, const Op& op) {
    typedef typename VecTraits<T>::reg_type reg_type;
    Vector<T> B (VecTraits<T>::stride,b);
    const reg_type* va = (const reg_type*) &a[0];
    const reg_type* vb = (const reg_type*) &B[0];
    reg_type* vc = (reg_type*) &c[0];
    size_t simd_n = std::floor(a.size()/VecTraits<T>::stride),
      start_r = simd_n*VecTraits<T>::stride;
    for (size_t i = 0; i < simd_n; ++i)
        vc[i] = op.packed(va[i], *vb);
    std::transform(a.begin()+start_r, a.end(), c.begin()+start_r, std::bind2nd(op,b));
}

#endif /* SRC_MATRIX_SIMDTRAITS_HPP_ */
