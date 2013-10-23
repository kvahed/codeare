/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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
		LINEAR, POLYNOMIAL, CSPLINE, CSPLINE_PERIODIC, AKIMA, AKIMA_PERIODIC
	};
}


template<class T> struct ITraits;

template<> struct ITraits<float> {
	typedef float T;
	static void ac (const container<T>& data, std::vector<container<double> >& y){
		size_t n = data.size();
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i)
			y[0][i] = data[i];
	}
	static void ac (const T* data, const size_t n, std::vector<container<double> >& y){
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i)
			y[0][i] = data[i];
	}
};
template<> struct ITraits<double> {
	typedef double T;
	static void ac (const container<T>& data, std::vector<container<double> >& y){
		y.push_back(data);
	}
	static void ac (const T* data, const size_t n, std::vector<container<double> >& y){
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i)
			y[0][i] = data[i];
	}
};
template<> struct ITraits<cxfl> {
	typedef cxfl T;
	static void ac (const container<T>& data, std::vector<container<double> >& y){
		size_t n = data.size();
		y.push_back(container<double>(n));
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i) {
			y[0][i] = real(data[i]); y[1][i] = imag(data[i]);
		}
	}
	static void ac (const T* data, const size_t n, std::vector<container<double> >& y){
		y.push_back(container<double>(n));
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i) {
			y[0][i] = real(data[i]); y[1][i] = imag(data[i]);
		}
	}
};
template<> struct ITraits<cxdb> {
	typedef cxdb T;
	static void ac (const container<T>& data, std::vector<container<double> >& y){
		y.push_back(real(data));
		y.push_back(imag(data));
	}
	static void ac (const T* data, const size_t n, std::vector<container<double> >& y){
		y.push_back(container<double>(n));
		y.push_back(container<double>(n));
		for (size_t i = 0; i < n; ++i) {
			y[0][i] = real(data[i]); y[1][i] = imag(data[i]);
		}
	}
};


/**
 * @brief    Evaluate polynomial
 */
template <class T>
class PolyVal {
	
public:

	
	/**
	 * @brief     Default constructor
	 */
	PolyVal () {};


	/**
	 * @brief     Copy constructor
	 */
	PolyVal (const PolyVal& pv) {
		*this = pv;
	}


	/**
	 * @brief     Assignement operator
	 */
	PolyVal& operator= (const PolyVal& pv) {
		if (this != &pv) {
			m_y = pv.m_y;
			m_acc = pv.m_acc;
			m_spline = pv.m_spline;
		}
		return *this;
	}

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
		
		assert (x.Size() == y.Size());

		if (!Initialise (x.Container(), y.Container(), intm))
			printf ("  PolyVal construction failed!\n");
		
	}

	
	/**
	 * @brief        Construct
	 * 
	 * @param  x     Vector of x (Matrix<double>)
	 * @param  y     Vector of y(x) (Matrix<double>)
	 * @param  intm  Interpolation Method (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 *
	 * @return       Success
	 */ 
	PolyVal (const container<double>& x, const container<T>& y, const INTERP::Method intm = INTERP::CSPLINE) {
		
		assert (x.size() == y.size());

		if (!Initialise (x, y, intm))
			printf ("  PolyVal construction failed!\n");
		
	}

	
	/**
	 * @brief        Construct
	 *
	 * @param  x     Vector of x (Matrix<double>)
	 * @param  y     Vector of y(x) (Matrix<double>)
	 * @param  intm  Interpolation Method (linear, polynomial, [periodic] cubic splice, [periodic] akima)
	 *
	 * @return       Success
	 */
	PolyVal (const Matrix<double>& x, const T* y, const INTERP::Method intm = INTERP::CSPLINE) {

		if (!Initialise (x.Container(), y, intm))
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
	PolyVal (const Matrix<double>& x, const container<T>& y, const INTERP::Method intm = INTERP::CSPLINE) {

		assert (x.Size() == y.size());

		if (!Initialise (x.Container(), &y[0], intm))
			printf ("  PolyVal construction failed!\n");

	}


	/**
	 * @brief    Clean up
	 */
	virtual ~PolyVal () {

		size_t i;

		for (i = 0; i < m_spline.size(); ++i)
				gsl_spline_free (m_spline[i]);
		for (i = 0; i < m_acc.size(); ++i)
				gsl_interp_accel_free (m_acc[i]);
		
	}


	/**
	 * @brief     Interpolate at point
	 *
	 * @param xx  Point
	 * @return    Value
	 */
	inline T 
	Lookup (const double xx) const; 


protected:
	

	inline bool InitGSL (const container<double>& x, const INTERP::Method intm) {

		size_t n = x.size();

		for (size_t i = 0; i < m_y.size(); ++i) {
			m_acc.push_back(gsl_interp_accel_alloc());
			m_spline.push_back(gsl_spline_alloc (IntMeth(intm), n));
			gsl_spline_init (m_spline[i], x.memory(), m_y[i].memory(),n);
		}

		return true;

	}

	bool Initialise (const container<double>& x, const container<T>& y, const INTERP::Method intm) {

		ITraits<T>::ac (y, m_y);
		return InitGSL (x, intm);

	}

	bool Initialise (const container<double>& x, const T* y, const INTERP::Method intm) {

		ITraits<T>::ac(y, x.size(), m_y);
		return InitGSL (x, intm);

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
		
	std::vector<container<double> > m_y;
	std::vector<gsl_interp_accel*> m_acc;    /**< @brief Accelerators for real / imaginary interpolations */
	std::vector<gsl_spline*>       m_spline; /**< @brief Spline finction for real / imaginary      */
	
};

template<> inline double 
PolyVal<double>::Lookup (const double xx) const {
	return gsl_spline_eval (m_spline[0], xx, m_acc[0]);
}

template<> inline float 
PolyVal<float>::Lookup (const double xx) const {
	return (float) gsl_spline_eval (m_spline[0], xx, m_acc[0]);
}

template<> inline cxfl 
PolyVal<cxfl>::Lookup (const double xx) const {
	return cxfl(gsl_spline_eval (m_spline[0], xx, m_acc[0]), gsl_spline_eval (m_spline[1], xx, m_acc[1]));
}

template<> inline cxdb 
PolyVal<cxdb>::Lookup (const double xx) const {
	return cxdb(gsl_spline_eval (m_spline[0], xx, m_acc[0]), gsl_spline_eval (m_spline[1], xx, m_acc[1]));
}

#endif /*__POLY_VAL_HPP__*/
