/*
 * NoiseTraits.hpp
 *
 *  Created on: Jun 15, 2013
 *      Author: kvahed
 */

#ifndef NOISETRAITS_HPP_
#define NOISETRAITS_HPP_

#include <gsl/gsl_rng.h>
#include <limits>

template<class T> class RandTraits;

/*
template<> class NoiseTraits<float> {

	const gsl_rng_type* grt;
	gsl_rng* r;
    double s;

	gsl_rng_env_setup();

    //switch (D){
    //   case normal:  grt = gsl_rng_default; break;
            //case gaussian: s = 1.0; grt = gsl_ran_gaussian; break;
    //}

	grt = gsl_rng_default;
	r = gsl_rng_alloc (grt);
	gsl_rng_set (r, getticks()); // reseed

	if      (typeid(T) == typeid(float) || typeid(T) == typeid(double))
		while (i--)
			res[i] = 2.0 * gsl_rng_uniform (r) - 1.0;
	else if (typeid(T) == typeid(cxdb))
		while (i--) {
			((double*) &res[i])[0] = 2.0 * gsl_rng_uniform (r) - 1.0;
			((double*) &res[i])[1] = 2.0 * gsl_rng_uniform (r) - 1.0;
		}
	else if (typeid(T) == typeid(cxfl))
		while (i--) {
			((float*) &res[i])[0] = 2.0 * gsl_rng_uniform (r) - 1.0;
			((float*) &res[i])[1] = 2.0 * gsl_rng_uniform (r) - 1.0;
		}
	else if (typeid(T) == typeid(short))
		while (i--)
			res[i] = 2 * gsl_rng_uniform_int (r, SHRT_MAX) - SHRT_MAX;
	else if (typeid(T) == typeid(long))
		while (i--)
			res[i] = 2 * gsl_rng_uniform_int (r, INT_MAX) - INT_MAX;

	gsl_rng_free (r);


};
*/
#endif /* NOISETRAITS_HPP_ */
