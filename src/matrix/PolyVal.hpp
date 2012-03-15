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

	/**
	 * @brief  1D interpolation method
	 */ 
	enum Method {

		LINEAR,           /**< @brief linear */
		POLYNOMIAL,       /**< @brief polynomial */
		CSPLINE,          /**< @brief cubic spline */
		CSPLINE_PERIODIC, /**< @brief cubic spline with peridicity constraint at edges */
		AKIMA,            /**< @brief akima */
		AKIMA_PERIODIC    /**< @brief akima with peridocity constraint at edges */
		
	};

}

/**
 * @brief    Evaluate polynomial (GSL interpolation)
 */
template <class T>
class PolyVal {
	
public:
	
	
	/**
	 * @brief        Default constructor
	 */
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
	PolyVal (const Matrix<double>& x, const Matrix<T>& y, const INTERP::Method intm) {
		
		m_tmp_i = false;
		m_tmp   = false;

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
	PolyVal (const Matrix<double>& x, const T* y, const INTERP::Method intm = INTERP::CSPLINE);


	/**
	 * @brief    Clean up
	 */
	virtual ~PolyVal () {

		if (m_initialised)
			gsl_spline_free (m_spline);
	   
		if (m_allocated)
			gsl_interp_accel_free (m_acc);

		if (m_btmp)
			free (m_tmp);
		
		if (m_btmp_i)
			free (m_tmp_i);
		
	}


	/**
	 * @brief     Interpolate at point
	 *
	 * @param xx  Point
	 * @return    Value
	 */
	T 
	Lookup        (const double& xx) const;


private:
	

	/**
	 * @brief        Initialise 
	 * 
	 * @param  x     Vector of x
	 * @param  intm  Interpolation Mmthod (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 * @param  n     Size of x/y(x)
	 *
	 * @return       Success
	 */ 
	inline bool Initialise (const double* x, const INTERP::Method& intm, const size_t& n) {

		bool cx = (typeid(T) == typeid(cxfl)) || (typeid(T) == typeid(cxdb));
		
		m_allocated     = false;
		m_initialised   = false;		
		m_allocated_i   = false;
		m_initialised_i = false;
		
		m_acc = gsl_interp_accel_alloc ();
		if (cx)
			m_acc_i = gsl_interp_accel_alloc ();
		
		m_spline = gsl_spline_alloc (InterpMethod(intm), n);
		if (cx) {
			m_spline_i    = gsl_spline_alloc (InterpMethod(intm), n);
			m_allocated_i = true;
		}
		m_allocated   = true;
		
		if        (typeid(T) == typeid(double) || typeid(T) == typeid(float)) {
			gsl_spline_init (m_spline, x, m_tmp, n);
		} else if (typeid(T) == typeid(cxfl)   || typeid(T) == typeid(cxfl) ) {
			gsl_spline_init (m_spline, x, m_tmp, n);
			gsl_spline_init (m_spline_i, x, m_tmp_i, n);
			m_initialised_i = true;
		}
		m_initialised = true;

		return true;
		
	} 
	

	const gsl_interp_type* 
	InterpMethod (const INTERP::Method& im) const {
		
		switch (im)
			{
			case INTERP::LINEAR:           return gsl_interp_linear;           break;
			case INTERP::POLYNOMIAL: 	   return gsl_interp_polynomial;      break;
			case INTERP::CSPLINE:          return gsl_interp_cspline;          break;
			case INTERP::CSPLINE_PERIODIC: return gsl_interp_cspline_periodic; break;
			case INTERP::AKIMA:            return gsl_interp_akima;            break;
			case INTERP::AKIMA_PERIODIC:   return gsl_interp_akima_periodic;   break;
			}
		
		return gsl_interp_cspline;
		
	}


	bool              m_initialised; /**< @brief Initialised */
	bool              m_initialised_i; /**< @brief Initialised */
	bool              m_allocated;   /**< @brief GSL workspace allocated */
	bool              m_allocated_i;   /**< @brief GSL workspace allocated */
	bool              m_btmp;         /**< @brief Allocated m_tmp */  
	bool              m_btmp_i;       /**< @brief Allocated m_tmp_i*/  
	
	gsl_interp_accel* m_acc;         /**< @brief Accelerator */
	gsl_interp_accel* m_acc_i;       /**< @brief Accelerator (for imaginary) */
	gsl_spline*       m_spline  ;    /**< @brief Spline      */
	gsl_spline*       m_spline_i;    /**< @brief Spline      (for imaginary) */

	double* m_tmp;                     /**< @brief real of complex y */ 
	double* m_tmp_i;                   /**< @brief imaginary of complex y */
	
};



template<>
PolyVal<cxdb>::PolyVal (const Matrix<double>& x, const cxdb* y, const INTERP::Method intm) {

	size_t nx = 0;
	
	m_tmp    = (double*) malloc (nx * sizeof (double));
	m_tmp_i  = (double*) malloc (nx * sizeof (double));
	
	m_btmp   = true;
	m_btmp_i = true;

	for (size_t i = 0; i < nx; i++) {
		m_tmp  [i] = y[i].real();
		m_tmp_i[i] = y[i].imag();
	}
	
	if (!Initialise (x.Data(), intm, nx))
		printf ("  PolyVal construction failed!\n");

}


template<>
PolyVal<cxfl>::PolyVal (const Matrix<double>& x, const cxfl* y, const INTERP::Method intm) {

	size_t nx = 0;
	
	m_tmp    = (double*) malloc (nx * sizeof (double));
	m_tmp_i  = (double*) malloc (nx * sizeof (double));
	
	m_btmp   = true;
	m_btmp_i = true;

	for (size_t i = 0; i < nx; i++) {
		m_tmp  [i] = y[i].real();
		m_tmp_i[i] = y[i].imag();
	}

	if (!Initialise (x.Data(), intm, nx))
		printf ("  PolyVal construction failed!\n");

}


template<>
PolyVal<float>::PolyVal (const Matrix<double>& x, const float* y, const INTERP::Method intm) {

	size_t nx = sizeof (x,0);

	m_tmp = (double*) malloc  (nx * sizeof (double));
	m_btmp = true;

	for (size_t i = 0; i < nx; i++) {
		m_tmp[i] = y[i];
	}

	if (!Initialise (x.Data(), intm, nx))
		printf ("  PolyVal construction failed!\n");
	
}


template<>
PolyVal<double>::PolyVal (const Matrix<double>& x, const double* y, const INTERP::Method intm) {
	
	size_t nx = sizeof (x,0);

	m_tmp = (double*) malloc  (nx * sizeof (double));
	m_btmp = true;

	memcpy (m_tmp, y, nx);

	if (!Initialise (x.Data(), intm, nx))
		printf ("  PolyVal construction failed!\n");
	
}


template<> inline double 
PolyVal<double>::Lookup (const double& xx) const {

	return gsl_spline_eval (m_spline, xx, m_acc);
	
}

template<> inline float 
PolyVal<float>::Lookup (const double& xx) const {

	return gsl_spline_eval (m_spline, xx, m_acc);
	
}


template<> inline cxfl
PolyVal<cxfl>::Lookup (const double& xx) const {

	return cxfl(gsl_spline_eval (m_spline, xx, m_acc), gsl_spline_eval (m_spline_i, xx, m_acc_i));
	
}

template<> inline cxdb
PolyVal<cxdb>::Lookup (const double& xx) const {

	return cxfl(gsl_spline_eval (m_spline, xx, m_acc), gsl_spline_eval (m_spline_i, xx, m_acc_i));
	
}

#endif /*__POLY_VAL_HPP__*/
