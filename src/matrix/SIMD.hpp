#ifndef __SIMD_HPP__
#define __SIMD_HPP__

#include "SIMDTraits.hpp"

#include <math.h>
#include <stdio.h>
#include <valarray>

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
	
#endif

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
	process (const std::valarray<T>& A, const std::valarray<T>& B, const Op& op, std::valarray<T>& C) {
		
		typedef SSETraits<T>                sse_type;
		typedef typename sse_type::Register reg_type;
		
		size_t   n  = A.size();
		size_t   ne = sse_type::ne;
		size_t   na = floor ((float)n / (float)ne);
		size_t   nr = n - na*ne;// sse_type::ns;
		
		load<MA>  ld;
		store<MA> st;
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
		for (size_t i = 0; i < nr; i++) {
			a = load<MA>::unaligned (pA++);
			b = load<MA>::unaligned (pB++);
			c = op.single (a, b);
			store<MA>::unaligned (pC++, c);
		}

	}; // namespace SSE


	
}

#endif // __SIMD_HPP__
