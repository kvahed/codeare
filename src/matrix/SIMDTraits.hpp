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
#include <immintrin.h>

template<class T> struct AVXTraits;
template<> struct AVXTraits<float> {
    typedef __m256 reg_type;
    static const int stride = 8;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_ps(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm256_mul_ps(a, b);}
};
template<> struct AVXTraits<double> {
    typedef __m256d reg_type;
    static const int stride = 4;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_pd(a, b);}
    inline static reg_type multiplies (const reg_type& a, const reg_type& b) {return _mm256_mul_pd(a, b);}
};
template<> struct AVXTraits<cxfl> {
    typedef __m256 reg_type;
    static const int stride = 4;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_ps(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_ps(a, b);}
    inline static reg_type multiplies (const reg_type &a, const reg_type &b) {
    	reg_type b_flip = _mm256_shuffle_ps(b,b,0xB1);   // Swap b.re and b.im
    	reg_type a_im   = _mm256_shuffle_ps(a,a,0xF5);   // Imag part of a in both
    	reg_type a_re   = _mm256_shuffle_ps(a,a,0xA0);   // Real part of a in both
    	reg_type aib    = _mm256_mul_ps(a_im, b_flip);   // (a.im*b.im, a.im*b.re)
#ifdef __FMA__      // FMA3
    	return  _mm256_fmaddsub_ps(a_re, b, aib);      // a_re * b +/- aib
#elif defined (__FMA4__)  // FMA4
    	return  _mm256_maddsub_ps (a_re, b, aib);      // a_re * b +/- aib
#else
    	reg_type arb    = _mm256_mul_ps(a_re, b);        // (a.re*b.re, a.re*b.im)
    	return          _mm256_addsub_ps(arb, aib);    // subtract/add
#endif  // FMA
    }
};
template<> struct AVXTraits<cxdb> {
    typedef __m256d reg_type;
    static const int stride = 2;
    inline static reg_type plus (const reg_type& a, const reg_type& b) {return _mm256_add_pd(a, b);}
    inline static reg_type minus (const reg_type& a, const reg_type& b) {return _mm256_sub_pd(a, b);}
    static inline reg_type multiplies (const reg_type &a, const reg_type &b) {
        reg_type b_flip = _mm256_shuffle_pd(b,b,5); // Swap b.re and b.im
        reg_type a_im = _mm256_shuffle_pd(a,a,0xF); // Imag part of a in both
        reg_type a_re = _mm256_shuffle_pd(a,a,0); // Real part of a in both
        reg_type aib = _mm256_mul_pd(a_im, b_flip); // (a.im*b.im, a.im*b.re)
#ifdef __FMA__ // FMA3
        return _mm256_fmaddsub_pd(a_re, b, aib); // a_re * b +/- aib
#elif defined (__FMA4__) // FMA4
        return _mm256_maddsub_pd (a_re, b, aib); // a_re * b +/- aib
#else
        reg_type arb = _mm256_mul_pd(a_re, b); // (a.re*b.re, a.re*b.im)
        return _mm256_addsub_pd(arb, aib); // subtract/add
#endif // FMA
    }

};

namespace codeare {
template<class T> class plus {
public: 
    typedef typename AVXTraits<T>::reg_type reg_type;
    inline static reg_type packed (const reg_type& a, const reg_type& b) {
        return AVXTraits<T>::plus(a, b);
    }
    inline T operator()( const T& x, const T& y ) const {
        return std::plus<T>()( x, y );
    }
};
template<class T> class minus {
public: 
    typedef typename AVXTraits<T>::reg_type reg_type;
    inline static reg_type packed (const reg_type& a, const reg_type& b) {
        return AVXTraits<T>::minus(a, b);
    }
    inline T operator()( const T& x, const T& y ) const {
        return std::minus<T>()( x, y );
    }
};
template<class T> class multiplies {
public: 
    typedef typename AVXTraits<T>::reg_type reg_type;
    inline static reg_type packed (const reg_type& a, const reg_type& b) {
        return AVXTraits<T>::multiplies(a, b);
    }
    inline T operator()( const T& x, const T& y ) const {
        return std::multiplies<T>()( x, y );
    }
};
}
template<class T, class Op> 
inline static void Vec (const Vector<T>& a, const Vector<T>& b, Vector<T>& c, const Op& op) {
    typedef typename TypeTraits<T>::RT RT;
    typedef typename AVXTraits<T>::reg_type reg_type;
    const reg_type* va = (const reg_type*) &a[0];
    const reg_type* vb = (const reg_type*) &b[0];
    reg_type* vc = (reg_type*) &c[0];
    int simd_n = std::floor(a.size()/AVXTraits<T>::stride);
    size_t start_r = simd_n*AVXTraits<T>::stride;
#pragma omp parallel for
    for (int i = 0; i < simd_n; ++i)
        vc[i] = op.packed(va[i], vb[i]);
    std::transform(a.begin()+start_r, a.end(), b.begin()+start_r, c.begin()+start_r, op);
}

#endif /* SRC_MATRIX_SIMDTRAITS_HPP_ */
