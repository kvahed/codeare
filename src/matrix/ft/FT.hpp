/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
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

#ifndef __FT_HPP__
#define __FT_HPP__

#include "CX.hpp"
#include "Params.hpp"
#include "Operator.hpp"

/**
 * @brief  Base class for single and double precision complex Fourier transforms
 */
template <class T>
class FT : public Operator<T> {

	typedef typename TypeTraits<T>::RT RT;

public:

	/**
	 * @brief    Default constructor
	 */
	FT () {};


	/**
	 * @brief     Contstruct with parameters
	 */
	FT (const Params& params) {
		m_params = params;
	}

	/**
	 * @brief    Default destructor
	 */
	virtual ~FT() {};

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
//	virtual Matrix<T> Trafo (const Matrix<T>& m) const { return Matrix<T>();}
	virtual Matrix<T> Trafo (const MatrixType<T>& m) const {return Matrix<T>();}
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
//	virtual Matrix<T> Adjoint (const Matrix<T>& m) const {return Matrix<T>();}
	virtual Matrix<T> Adjoint (const MatrixType<T>& m) const {return Matrix<T>();}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> operator* (const Matrix<T>& m) const {
		return Trafo(m);
	}
	virtual Matrix<T> operator* (const MatrixType<T>& m) const {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> operator->* (const Matrix<T>& m) const {
		return Adjoint (m);
	}
	virtual Matrix<T> operator->* (const MatrixType<T>& m) const {
		return Adjoint (m);
	}

    virtual FT<T>* getFT() { return this; }

	/**
	 * @brief      Assign k-space trajectory
	 *
	 * @param  k   K-space trajectory
	 */
	virtual void KSpace (const Matrix<RT>& k) {}

	/**
	 * @brief      Assign k-space weigths (jacobian of k in t)
	 *
	 * @param  w   Weights
	 */
	virtual void Weights (const Matrix<RT>& w) {}

	/**
	 * @brief      Assign k-space weigths (jacobian of k in t)
	 *
	 * @param  w   Weights
	 */
	virtual void Mask (const Matrix<RT>& m) {}


protected:
	
	Params m_params;

};

#endif
