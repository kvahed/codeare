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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <limits>

template<class T> inline static Matrix<T>
normal (Matrix<T>&) {}


enum distribution { uniform, gaussian, rayleigh, landau };
enum algorithm { def, mt19937, ranlxs0, ranlxs1, ranlxs2, ranlxd1, ranlxd2,
	             ranlux, ranlux389, cmrg, mrg, taus, taus2, gfsr4};

inline static const gsl_rng_type*
gen_type (const algorithm alg) {

	const gsl_rng_type* T;

	switch (alg) {
		case 0: T=gsl_rng_default;   break;
		case 1: T=gsl_rng_mt19937;   break;
		case 2: T=gsl_rng_ranlxs0;   break;
		case 3: T=gsl_rng_ranlxs1;   break;
		case 4: T=gsl_rng_ranlxs2;   break;
		case 5: T=gsl_rng_ranlxd1;   break;
		case 6: T=gsl_rng_ranlxd2;   break;
		case 7: T=gsl_rng_ranlux;    break;
		case 8: T=gsl_rng_ranlux389; break;
		case 9: T=gsl_rng_cmrg;      break;
		case 10: T=gsl_rng_mrg;      break;
		case 11: T=gsl_rng_taus;     break;
		case 12: T=gsl_rng_taus2;    break;
		case 13: T=gsl_rng_gfsr4;    break;
		default: T=gsl_rng_default;  break;
	}

	return T;

}

template <class T, distribution D>
struct generator;

template<> template<>
struct generator<float,uniform> {
	typedef float T;
	inline static T generate (gsl_rng* r) {
		return 2.0 * gsl_rng_uniform (r) - 1.0;
	}
};
template<> template<>
struct generator<double,uniform> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return 2.0 * gsl_rng_uniform (r) - 1.0;
	}
};
template<> template<>
struct generator<cxfl,uniform> {
	typedef cxfl T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_rng_uniform(r) - 1.0, 2.0 * gsl_rng_uniform(r) - 1.0);
	}
};
template<> template<>
struct generator<cxdb,uniform> {
	typedef cxdb T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_rng_uniform(r) - 1.0, 2.0 * gsl_rng_uniform(r) - 1.0);
	}
};
template<> template<>
struct generator<short,uniform> {
	typedef short T;
	inline static T generate (gsl_rng* r) {
		return 2 * gsl_rng_uniform_int (r, SHRT_MAX) - SHRT_MAX;
	}
};
template<> template<>
struct generator<long,uniform> {
	typedef long T;
	inline static T generate (gsl_rng* r) {
		return 2 * gsl_rng_uniform_int (r, INT_MAX) - INT_MAX;
	}
};
template<> template<>
struct generator<unsigned short,uniform> {
	typedef unsigned short T;
	inline static T generate (gsl_rng* r) {
		return gsl_rng_uniform_int (r, SHRT_MAX);
	}
};
template<> template<>
struct generator<size_t,uniform> {
	typedef size_t T;
	inline static T generate (gsl_rng* r) {
		return gsl_rng_uniform_int (r, INT_MAX);
	}
};

template<> template<>
struct generator<float,gaussian> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_ugaussian (r) - 1.0);
	}
};
template<> template<>
struct generator<double,gaussian> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_ugaussian (r) - 1.0);
	}
};
template<> template<>
struct generator<cxfl,gaussian> {
	typedef cxfl T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_ugaussian (r) - 1.0, 2.0 * gsl_ran_ugaussian (r) - 1.0);
	}
};
template<> template<>
struct generator<cxdb,gaussian> {
	typedef cxdb T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_ugaussian (r) - 1.0, 2.0 * gsl_ran_ugaussian (r) - 1.0);
	}
};
template<> template<>
struct generator<float,rayleigh> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_rayleigh (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<double,rayleigh> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_rayleigh (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<cxfl,rayleigh> {
	typedef cxfl T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_rayleigh (r,1.0) - 1.0, 2.0 * gsl_ran_rayleigh (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<cxdb,rayleigh> {
	typedef cxdb T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_rayleigh (r,1.0) - 1.0, 2.0 * gsl_ran_rayleigh (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<float,landau> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_landau (r) - 1.0);
	}
};
template<> template<>
struct generator<double,landau> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_landau (r) - 1.0);
	}
};
template<> template<>
struct generator<cxfl,landau> {
	typedef cxfl T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_landau (r) - 1.0, 2.0 * gsl_ran_landau (r) - 1.0);
	}
};
template<> template<>
struct generator<cxdb,landau> {
	typedef cxdb T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_landau (r) - 1.0, 2.0 * gsl_ran_landau (r) - 1.0);
	}
};

template<class T, distribution D>
inline static void rand_pop (Matrix<T>& M, algorithm a = def) {

	const gsl_rng_type* grt = gen_type(a);
	gsl_rng_env_setup();
	gsl_rng* r = gsl_rng_alloc (grt);

	gsl_rng_set (r, getticks()); // reseed

	for (size_t i = 0; i < numel(M); ++i)
		M[i] = generator<T,D>::generate(r);

	gsl_rng_free (r);

}


#include <boost/random.hpp> 




template<class T, class G = boost::mt19937, class D = boost::normal_distribution<typename TypeTraits<T>::RT > >
class Random {

    typedef boost::variate_generator<G&,D> vg_type;
    
public:
    
    inline static void Populate (Matrix<T>& M, const T mean, const T sigma) {
        static G rng (static_cast<unsigned> (time(0)));
        vg_type generator (rng, D(mean,sigma));
        std::generate (M.Container().begin(), M.Container.end(), generator);
    }
    
};


template<class T> inline std::complex<T>
AddNoise (std::complex<T> &val, double &noiseest, const double relativenoiselevel,
          const double absolutenoiselevel, boost::lagged_fibonacci607 &generator) {
    
    typedef boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> > vg_type;
    typedef boost::normal_distribution<> d_type;

    //determine the noise level for the real and imaginary parts
    double realnoiselevel = std::max (std::abs(val.real() * relativenoiselevel), absolutenoiselevel);
    double imagnoiselevel = std::max (std::abs(val.imag() * relativenoiselevel), absolutenoiselevel);
    
    //draw a sample from a normal distribution
    return std::complex<T>(vg_type (generator, dtype(val.real(), realnoiselevel))(),
                           vg_type (generator, dtype(val.imag(), imagnoiselevel))());
    
}

#endif /* NOISETRAITS_HPP_ */
