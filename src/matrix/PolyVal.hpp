/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#ifndef __POLY_VAL_HPP__
#define __POLY_VAL_HPP__

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Algos.hpp"

namespace INTERP {

	enum Method {

		LINEAR,
		POLYNOMIAL,
		CSPLINE,
		CSPLINE_PERIODIC,
		AKIMA,
		AKIMA_PERIODIC
		
	};

}

/**
 * @brief    Evaluate polynomial
 */
template <class T>
class PolyVal {
	
public:

	
	PolyVal () {

		Reset();

	};


	/**
	 * @brief        Construct
	 * 
	 * @param  x     Vector of x (Matrix<double>)
	 * @param  y     Vector of y(x) (Matrix<double>)
	 * @param  intm  Interpolation Method (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 *
	 * @return       Success
	 */ 
	PolyVal (const Matrix<double>& x, const Matrix<T>& y, const INTERP::Method intm = INTERP::CSPLINE) {

		Reset();
		
		if (!Initialise (x.Data(), y.Data(), intm, size(x,0)))
			printf ("  PolyVal construction failed!\n");
		
	}
	
	
	/**
	 * @brief        Construct
	 * 
	 * @param  x     Vector of x (Matrix<double>)
	 * @param  y     Vector of y(x) (double*)
	 * @param  intm  Interpolation Method (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 *
	 * @return       Success
	 */ 
	PolyVal (const Matrix<double>& x, const T* y, const INTERP::Method intm = INTERP::CSPLINE) {
	
		Reset ();
	
		if (!Initialise (x.Data(), y, intm, size(x,0)))
			printf ("  PolyVal construction failed!\n");

	}


	inline void Reset () {

		m_alloc[0] = false;
		m_alloc[1] = false;
		m_init [0] = false;
		m_init [1] = false;

	}


	/**
	 * @brief    Clean up
	 */
	virtual ~PolyVal () {

		for (int i = 0; i < 2; i++) {
			
			if (m_init[i])
				gsl_spline_free (m_spline[i]);
			
			if (m_alloc[i]) {
				gsl_interp_accel_free (m_acc[i]);
				free (m_y[i]);
			}

		}
		
	}


	/**
	 * @brief     Interpolate at point
	 *
	 * @param xx  Point
	 * @return    Value
	 */
	T Lookup (const double& xx) const; 


private:
	

	/**
	 * @brief        Initialise 
	 * 
	 * @param  x     Vector of x
	 * @param  y     Vector of y(x)
	 * @param  intm  Interpolation Mmthod (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 * @param  n     Size of x/y(x)
	 *
	 * @return       Success
	 */ 
	bool Initialise (const double* x, const T* y, const INTERP::Method intm, const size_t& n);


	inline bool InitGSL (const double* x, const INTERP::Method intm, const size_t& n, const int& p) {

		m_acc   [p] = gsl_interp_accel_alloc ();
		m_spline[p] = gsl_spline_alloc (IntMeth(intm), n);
		m_alloc [p] = true;
		
		gsl_spline_init (m_spline[p], x, m_y[p], n);
		m_init  [p] = true;
		
		return (m_init [p] && m_alloc [p]);

	}


	inline const gsl_interp_type* 
	IntMeth (const INTERP::Method& intm) const {

		switch (intm)
			{
			case INTERP::LINEAR:           return gsl_interp_linear; break;
			case INTERP::POLYNOMIAL: 	   return gsl_interp_polynomial; break;
			case INTERP::CSPLINE:          return gsl_interp_cspline; break;
			case INTERP::CSPLINE_PERIODIC: return gsl_interp_cspline_periodic; break;
			case INTERP::AKIMA:            return gsl_interp_akima; break;
			case INTERP::AKIMA_PERIODIC:   return gsl_interp_akima_periodic; break;
			}

		return gsl_interp_cspline;

	}
		
	double*           m_y     [2]; /**< @brief All types but doubles need conversion */

	
	
	bool              m_init  [2]; /**< @brief Real / imaginary interpolation initialised */
	bool              m_alloc [2]; /**< @brief Real / imaginary interpolation allocated */

	gsl_interp_accel* m_acc   [2];    /**< @brief Accelerators for real / imaginary interpolations */
	gsl_spline*       m_spline[2]; /**< @brief Spline finction for real / imaginary      */
	
};


template<> inline bool 
PolyVal<double>::Initialise (const double* x, const double* y, const INTERP::Method intm, const size_t& n) {
	
	size_t nn = n * sizeof (double);
	
	m_y [0] = (double*) malloc (nn);
	memcpy (m_y[0], y, nn);
	
	return InitGSL (x, intm, n, 0);
	
} 


template<> inline bool 
PolyVal<float>::Initialise (const double* x, const float* y, const INTERP::Method intm, const size_t& n) {
	
	size_t nn = n * sizeof (double);
	
	m_y [0] = (double*) malloc (nn);

	for (size_t i = 0; i < n; i++)
		m_y[0][i] = y[i];
	
	return InitGSL (x, intm, n, 0);
	
} 


template<> inline bool 
PolyVal<cxdb>::Initialise (const double* x, const cxdb* y, const INTERP::Method intm, const size_t& n) {
	
	size_t nn = n * sizeof (double);
	
	m_y [0] = (double*) malloc (nn);
	m_y [1] = (double*) malloc (nn);

	for (size_t i = 0; i < n; i++) {
		m_y[0][i] = y[i].real();
		m_y[1][i] = y[i].imag();
	}
	
	bool ri = InitGSL (x, intm, n, 0);
	bool ii = InitGSL (x, intm, n, 1);
	
	return (ri && ii);
	
} 


template<> inline bool 
PolyVal<cxfl>::Initialise (const double* x, const cxfl* y, const INTERP::Method intm, const size_t& n) {
	
	size_t nn = n * sizeof (double);
	
	m_y [0] = (double*) malloc (nn);
	m_y [1] = (double*) malloc (nn);

	for (size_t i = 0; i < n; i++) {
		m_y[0][i] = y[i].real();
		m_y[1][i] = y[i].imag();
	}
	
	bool ri = InitGSL (x, intm, n, 0);
	bool ii = InitGSL (x, intm, n, 1);
	
	return (ri && ii);
	
} 



template<> inline double 
PolyVal<double>::Lookup (const double& xx) const {
	
	return gsl_spline_eval (m_spline[0], xx, m_acc[0]);
	
}

template<> inline float 
PolyVal<float>::Lookup (const double& xx) const {
	
	return (float) gsl_spline_eval (m_spline[0], xx, m_acc[0]);
	
}

template<> inline cxfl 
PolyVal<cxfl>::Lookup (const double& xx) const {
	
	return cxfl(gsl_spline_eval (m_spline[0], xx, m_acc[0]), gsl_spline_eval (m_spline[1], xx, m_acc[1]));
	
}

template<> inline cxdb 
PolyVal<cxdb>::Lookup (const double& xx) const {
	
	return cxdb(gsl_spline_eval (m_spline[0], xx, m_acc[0]), gsl_spline_eval (m_spline[1], xx, m_acc[1]));
	
}

#endif /*__POLY_VAL_HPP__*/
