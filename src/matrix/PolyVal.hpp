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

	
	PolyVal () : m_initialised (false) {};


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
		
		if (!Initialise (x.Data(), y, intm, size(x,0)))
			printf ("  PolyVal construction failed!\n");

	}


	/**
	 * @brief    Clean up
	 */
	virtual ~PolyVal () {

		if (m_initialised)
			gsl_spline_free (m_spline);
	   
		if (m_allocated)
			gsl_interp_accel_free (m_acc);

		
	}


	/**
	 * @brief     Interpolate at point
	 *
	 * @param xx  Point
	 * @return    Value
	 */
	inline T Lookup (const double& xx) const {

		return gsl_spline_eval (m_spline, xx, m_acc);

	}


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
	inline bool Initialise (const double* x, const T* y, const INTERP::Method intm, const size_t n) {

		m_allocated   = false;
		m_initialised = false;
		
		m_acc = gsl_interp_accel_alloc ();
		
		switch (intm)
			{
			case INTERP::LINEAR:           m_spline = gsl_spline_alloc (gsl_interp_linear,            n); break;
			case INTERP::POLYNOMIAL: 	   m_spline = gsl_spline_alloc (gsl_interp_polynomial,        n); break;
			case INTERP::CSPLINE:          m_spline = gsl_spline_alloc (gsl_interp_cspline,           n); break;
			case INTERP::CSPLINE_PERIODIC: m_spline = gsl_spline_alloc (gsl_interp_cspline_periodic,  n); break;
			case INTERP::AKIMA:            m_spline = gsl_spline_alloc (gsl_interp_akima,             n); break;
			case INTERP::AKIMA_PERIODIC:   m_spline = gsl_spline_alloc (gsl_interp_akima_periodic,    n); break;
			default:                       printf ("GSL: INTERP ERROR - Invalid method.\n"); return false; break;
			}
		
		m_allocated   = true;

		gsl_spline_init (m_spline, x, y, n);
		m_initialised = true;

		return true;

	} 


	bool              m_initialised;
	bool              m_allocated;

	gsl_interp_accel* m_acc;      /**< @brief Accelerator */
	gsl_interp_accel* m_acc_i;    /**< @brief Accelerator (for imaginary) */
	gsl_spline*       m_spline  ; /**< @brief Spline      */
	gsl_spline*       m_spline_i; /**< @brief Spline      (for imaginary) */
	
};

#endif /*__POLY_VAL_HPP__*/
