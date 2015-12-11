/*
 * NoiseTraits.hpp
 *
 *  Created on: Jun 15, 2013
 *      Author: kvahed
 */

#ifndef NOISETRAITS_HPP_
#define NOISETRAITS_HPP_

#include "TypeTraits.hpp"
#include "Matrix.hpp"
#include <time.h>
#include <limits>

#include <limits>
#include <iostream>
#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>

template <class T>
struct RandTraits;

template <>
struct RandTraits<float> {
    typedef float T;
    typedef TypeTraits<T>::RT RT;
    typedef boost::normal_distribution<RT> normal;
    typedef boost::random::uniform_real_distribution<RT> uniform;
    inline static RT stdmin () {return -1.0;}
    inline static RT stdmax () {return +1.0;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) {
    	std::cout << "+++" << rta.Size() << std::endl;
		for (size_t i = 0; i < rta.Size(); ++i)
			rta[i] = generator();
	}
};
template <>
struct RandTraits<double> {
    typedef double T;
    typedef  TypeTraits<T>::RT RT;
    typedef boost::normal_distribution<RT> normal;
    typedef boost::random::uniform_real_distribution<RT> uniform;
    inline static RT stdmin () {return -1.0;}
    inline static RT stdmax () {return +1.0;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) {
		for (size_t i = 0; i < rta.Size(); ++i)
			rta[i] = generator();
	}
};
template <>
struct RandTraits<cxfl> {
    typedef cxfl T;
    typedef  TypeTraits<T>::RT RT;
    typedef boost::random::normal_distribution<RT> normal;
    typedef boost::random::uniform_real_distribution<RT> uniform;
    inline static RT stdmin () {return -1.0;}
    inline static RT stdmax () {return +1.0;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) { 
        for (size_t i = 0; i < rta.Size(); ++i)
          rta[i] = T(generator(),generator());
    }
};
template <>
struct RandTraits<cxdb> {
    typedef cxdb T;
    typedef  TypeTraits<T>::RT RT;
    typedef boost::normal_distribution<RT> normal;
    typedef boost::random::uniform_real_distribution<RT> uniform;
    inline static RT stdmin () {return -1.0;}
    inline static RT stdmax () {return +1.0;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) { 
        for (size_t i = 0; i < rta.Size(); ++i)
          rta[i] = T(generator(),generator());
    }
};
template <>
struct RandTraits<short> {
    typedef short T;
    typedef  TypeTraits<T>::RT RT;
    typedef boost::random::uniform_int_distribution<RT> uniform;
    typedef boost::random::uniform_int_distribution<RT> normal;
    inline static RT stdmin () {return -SHRT_MAX;}
    inline static RT stdmax () {return +SHRT_MAX;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) {
        std::generate (rta.Begin(), rta.End(), generator);
    }
};
template <>
struct RandTraits<long> {
    typedef long T;
    typedef  TypeTraits<T>::RT RT;
    typedef boost::random::uniform_int_distribution<RT> uniform;
    typedef boost::random::uniform_int_distribution<RT> normal;
    inline static RT stdmin () {return -LONG_MAX;}
    inline static RT stdmax () {return +LONG_MAX;}
    template <class G>
    inline static void Populate (Matrix<T>& rta, G& generator) {
        std::generate (rta.Begin(), rta.End(), generator);
    }
};

template<class T, class RNG = boost::mt19937>
class Random {

    typedef typename TypeTraits<T>::RT RT;
    
public:
    
    inline static void Normal (Matrix<T>& vt, const RT mean = 0.0, const RT sigma = 1.0) {
        RNG rng (static_cast<unsigned> ((unsigned)clock()));
        boost::variate_generator<RNG&,typename RandTraits<T>::normal>
            generator (rng, typename RandTraits<T>::normal(mean,sigma));
        RandTraits<T>::Populate(vt, generator);
    }
    
    inline static void Uniform (Matrix<T>& vt, const RT min = RandTraits<T>::stdmin(),
                                const RT max = RandTraits<T>::stdmax()) {
        RNG rng (static_cast<unsigned> ((unsigned)clock()));
        boost::variate_generator<RNG&,typename RandTraits<T>::uniform>
            generator (rng, typename RandTraits<T>::uniform(min,max));
        RandTraits<T>::Populate(vt, generator);
    }

};

#endif /* NOISETRAITS_HPP_ */
