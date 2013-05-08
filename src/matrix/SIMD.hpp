#ifndef __SIMD_HPP__
#define __SIMD_HPP__

#include "config.h"
#include "SIMDTraits.hpp"
#include "AlignmentAllocator.hpp"

#include <vector>


#if defined USE_VALARRAY
    #include <valarray>
    #define VECTOR_TYPE(A) std::valarray<A>
    #define VECTOR_CONSTR(A,B) std::valarray<A>(B)
#else
    #include "AlignmentAllocator.hpp"
    #if defined __AVX__
        #define ALIGNEMENT 32
    #elif defined __SSE2__
        #define ALIGNEMENT 16
    #endif
    #define VECTOR_TYPE(A) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >
    #define VECTOR_CONSTR(A,B) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >(B)
    #define VECTOR_CONSTR_VAL(A,B,C) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >(B,C)
#endif


#include <math.h>
#include <stdio.h>

enum alignment {
	A,
	U
};

namespace SSE {

#ifdef HAVE_SSE

	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_load
	 */
	template<class T>
	class load {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		aligned (const T* p) { 
			return sse_type::loada (p); 
		}

		inline static reg_type 
		unaligned (const T* p) { 
			return sse_type::loadu (p); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_store
	 */
	template<class T>
	class store {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static void
		aligned (T* p, const reg_type& a) { 
			sse_type::stora (p, a); 
		}

		inline static void
	    unaligned (T* p, const reg_type& a) { 
			sse_type::storu (p, a); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_add
	 */
	template<class T>
	class add  {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::addp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::adds (a, b); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_sub
	 */
	template<class T>
	class sub {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::subp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::subs (a, b); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_mul
	 */
	template<class T>
	class mul {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::mulp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::muls (a, b); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_sqrt
	 */
	template<class T>
	class div {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::divp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::divs (a, b); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_sqrt
	 */
	template<class T>
	class sqrt {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a) { 
			return sse_type::sqrtp (a); 
		}

		inline static reg_type 
		single (const reg_type& a) { 
			return sse_type::sqrts (a); 
		}

	};
	
	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_min
	 */
	template<class T>
	class min {

	public:

		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::minp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::mins (a, b); 
		}

	};

	/**
	 * @brief  Templated interface to SSE/SSE2 _mm_max
	 */
	template<class T>
	class max {

	public:
		
		typedef typename SSETraits<T>::Register reg_type;
		typedef SSETraits<T> sse_type;

		inline static reg_type 
		packed (const reg_type& a, const reg_type& b) { 
			return sse_type::maxp (a, b); 
		}

		inline static reg_type 
		single (const reg_type& a, const reg_type& b) { 
			return sse_type::maxs (a, b); 
		}

	};
	

	/**
	 * @brief     Process SSE operation Op such that C = op (A, B);
	 *
	 * @param  A  Vector A
	 * @param  B  Vector B
	 * @param  N  Length of A, B and C
	 * @param  op Operator class
	 * @param  C  Vector C
	 */
	template<class MA, class T, class Op> inline static void
	binary (const VECTOR_TYPE(T)& A, const VECTOR_TYPE(T)& B, const Op& op, VECTOR_TYPE(T)& C) {
		
		typedef SSETraits<T>                sse_type;
		typedef typename sse_type::Register reg_type;
		
		size_t   n  = A.size();
		size_t   ne = sse_type::ne;
		size_t   na = floor ((float)n / (float)ne);
        size_t   nr = n % ne;
		
		reg_type a, b, c;
		
		const T* pA = &A[0];
		const T* pB = &B[0];
		T*       pC = &C[0];

		// aligned 

		for (size_t i = 0; i < na; i++) {
			a = load<MA>::aligned (pA);
			b = load<MA>::aligned (pB);
			c = op.packed (a, b);
			store<MA>::aligned (pC, c);
			pA += ne; pB += ne; pC += ne;
		}

		// rest
		if (nr) {
			a = load<MA>::unaligned (pA);
			b = load<MA>::unaligned (pB);
			c = op.single (a, b);
			store<MA>::unaligned (pC, c);
		}

	} // namespace SSE

	template<class MA, class T> inline static void
	assign (const VECTOR_TYPE(T)& A, VECTOR_TYPE(T)& B) {

		typedef SSETraits<T>                sse_type;
		typedef typename sse_type::Register reg_type;

		size_t   n  = A.size();
		size_t   ne = sse_type::ne;
		size_t   na = floor ((float)n / (float)ne);
        size_t   nr = n % ne;

		const T* pA = &A[0];
		      T* pB = &B[0];

		// aligned

		for (size_t i = 0; i < na; i++) {
			store<MA>::aligned (pB, load<MA>::aligned (pA));
			pA += ne; pB += ne;
		}

		// rest
		if (nr)
			store<MA>::unaligned (pB, load<MA>::unaligned (pA));

	}

	#endif

	
}

#endif // __SIMD_HPP__
