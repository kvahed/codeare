/*
 * NoiseTraits.hpp
 *
 *  Created on: Jun 15, 2013
 *      Author: kvahed
 */

#ifndef NOISETRAITS_HPP_
#define NOISETRAITS_HPP_

#include "Matrix.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <limits>

template<class T> inline static Matrix<T>
normal (Matrix<T>&) {}


enum distribution { uniform, gaussian };
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
		return (2.0 * gsl_ran_gaussian (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<double,gaussian> {
	typedef double T;
	inline static T generate (gsl_rng* r) {
		return (2.0 * gsl_ran_gaussian (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<cxfl,gaussian> {
	typedef cxfl T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_gaussian (r,1.0) - 1.0, 2.0 * gsl_ran_gaussian (r,1.0) - 1.0);
	}
};
template<> template<>
struct generator<cxdb,gaussian> {
	typedef cxdb T;
	inline static T generate (gsl_rng* r) {
		return T(2.0 * gsl_ran_gaussian (r,1.0) - 1.0, 2.0 * gsl_ran_gaussian (r,1.0) - 1.0);
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



#endif /* NOISETRAITS_HPP_ */
