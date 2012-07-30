#include <emmintrin.h>

namespace SSE {


	template<class T> 
	struct SSETraits {};


	template<>
	struct SSETraits<double> {
		
		typedef __m128d Register;         /**< @brief register type */
		static const unsigned int ne = 2; /**< @brief # of processed elements */
		
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
			return _mm_mul_pd(a, b);
		}
	
		/**
		 * @brief     SSE2 single multiplication
		 */
		static inline Register 
		muls (const Register &a, const Register &b) {
			return _mm_mul_sd(a, b);
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
		
	}; // SSETraits<double>
	

	template<>
	struct SSETraits<float> {
		
		typedef __m128 Register;         /**< @brief register type */
		static const unsigned int ne = 4; /**< @brief # of processed elements */
		
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

	
} // namespace SSE 
