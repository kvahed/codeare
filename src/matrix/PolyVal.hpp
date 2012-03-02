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

class PolyVal {
	
public:

	
	PolyVal () : m_initialised (false) {};

	PolyVal (Matrix<double>& x, Matrix<double>& y, const INTERP::Method intm = INTERP::CSPLINE) : m_initialised (false) {
		
		size_t nx;

		m_acc = gsl_interp_accel_alloc ();
		nx    = size(x,0);

		if      (intm == INTERP::LINEAR)
			m_spline = gsl_spline_alloc (gsl_interp_linear, nx);
		else if (intm == INTERP::POLYNOMIAL)
			m_spline = gsl_spline_alloc (gsl_interp_polynomial, nx);
		else if (intm == INTERP::CSPLINE)
			m_spline = gsl_spline_alloc (gsl_interp_cspline, nx);
		else if (intm == INTERP::CSPLINE_PERIODIC)
			m_spline = gsl_spline_alloc (gsl_interp_cspline_periodic, nx);
		else if (intm == INTERP::AKIMA)
			m_spline = gsl_spline_alloc (gsl_interp_akima, nx);
		else if (intm == INTERP::AKIMA_PERIODIC)
			m_spline = gsl_spline_alloc (gsl_interp_akima_periodic, nx);
	
		gsl_spline_init (m_spline, &x[0], &y[0], nx);

		m_initialised = true;
		
	}
	
	
	PolyVal (Matrix<double>& x, double* y, const INTERP::Method intm = INTERP::CSPLINE)  : m_initialised (false) {
		
		size_t nx;

		m_acc = gsl_interp_accel_alloc ();
		nx    = size(x,0);

		if      (intm == INTERP::LINEAR)
			m_spline = gsl_spline_alloc (gsl_interp_linear, nx);
		else if (intm == INTERP::POLYNOMIAL)
			m_spline = gsl_spline_alloc (gsl_interp_polynomial, nx);
		else if (intm == INTERP::CSPLINE)
			m_spline = gsl_spline_alloc (gsl_interp_cspline, nx);
		else if (intm == INTERP::CSPLINE_PERIODIC)
			m_spline = gsl_spline_alloc (gsl_interp_cspline_periodic, nx);
		else if (intm == INTERP::AKIMA)
			m_spline = gsl_spline_alloc (gsl_interp_akima, nx);
		else if (intm == INTERP::AKIMA_PERIODIC)
			m_spline = gsl_spline_alloc (gsl_interp_akima_periodic, nx);
	
		gsl_spline_init (m_spline, &x[0], y, nx);
		
		m_initialised = true;

	}


	virtual ~PolyVal () {

		if (m_initialised) {
			gsl_spline_free (m_spline);
			gsl_interp_accel_free (m_acc);
		}
		
	}


	inline double Lookup (const double& xx) const {

		return gsl_spline_eval (m_spline, xx, m_acc);

	}


private:
	
	bool              m_initialised;
	gsl_interp_accel* m_acc;    /**< @brief Accelerator */
	gsl_spline*       m_spline; /**< @brief Spline      */
	
};

#endif /*__POLY_VAL_HPP__*/
