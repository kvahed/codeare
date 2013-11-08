/*
 * SIMDTraits.hpp
 *
 *  Created on: May 3, 2013
 *      Author: kvahed
 */

#ifndef __SIMDTRAITS_HPP__
#define __SIMDTRAITS_HPP__

#include <complex>
#include <xmmintrin.h>


namespace SSE {


	#if defined __AVX__
		#include "AVXTraits.hpp"
	#elif defined __SSE2__
		#include "SSE2Traits.hpp"
	#endif

	//	#if defined __AVX2__
	//		#include "AVX2Traits.hpp"
//#elif defined __SSE4_2__
		#include "SSE42Traits.hpp"
//	#endif

} // namespace SSE


#endif /* __SIMDTRAITS_HPP__ */
