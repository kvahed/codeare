/*
 * SSE42Traits.hpp
 *
 *  Created on: May 3, 2013
 *      Author: kvahed
 */

#ifndef SSE42TRAITS_HPP_
#define SSE42TRAITS_HPP_

#include "SIMDTraits.hpp"

#include <smmintrin.h>
#include <nmmintrin.h>

template<>
struct SSETraits<size_t> {

    typedef __m128i Register;         /**< @brief register type */
    typedef size_t Type;
    static const unsigned int ne = 2; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register
    loada (const Type* p) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const Type* p) {
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register
    adds (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register
    subp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register
    subs (const Register &a, const Register &b) {
        return a;
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register
    mulp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register
    muls (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register
    divp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single division
     */
    static inline Register
    divs (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register
    sqrtp (const Register &a) {
    	return a;
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register
    sqrts (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    minp (const Register &a, const Register &b) {
    	return a;
   }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    mins (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    maxp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    maxs (const Register &a, const Register &b) {
    	return a;
   }
};

template<>
struct SSETraits<short> {

    typedef __m128i Register;         /**< @brief register type */
    typedef short Type;
    static const unsigned int ne = 2; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register
    loada (const Type* p) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const Type* p) {
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register
    adds (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register
    subp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register
    subs (const Register &a, const Register &b) {
        return a;
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register
    mulp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register
    muls (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register
    divp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single division
     */
    static inline Register
    divs (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register
    sqrtp (const Register &a) {
    	return a;
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register
    sqrts (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    minp (const Register &a, const Register &b) {
    	return a;
   }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    mins (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    maxp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    maxs (const Register &a, const Register &b) {
    	return a;
   }

};


template<>
struct SSETraits<long> {

    typedef __m128i Register;         /**< @brief register type */
    typedef long Type;
    static const unsigned int ne = 2; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register
    loada (const Type* p) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const Type* p) {
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (Type* p, Register a) {
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (Type* p, Register a) {
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register
    adds (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register
    subp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register
    subs (const Register &a, const Register &b) {
        return a;
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register
    mulp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register
    muls (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register
    divp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single division
     */
    static inline Register
    divs (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register
    sqrtp (const Register &a) {
    	return a;
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register
    sqrts (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    minp (const Register &a, const Register &b) {
    	return a;
   }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    mins (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register
    maxp (const Register &a, const Register &b) {
    	return a;
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register
    maxs (const Register &a, const Register &b) {
    	return a;
   }
};



#endif /* SSE42TRAITS_HPP_ */
