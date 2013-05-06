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


# ifndef __DWT_HPP__
# define __DWT_HPP__


/**
 * @brief  Supported wavelet families
 */
enum wlfamily {
	
	ID = -1,                  /**< Identity transform*/
	WL_DAUBECHIES,
	WL_DAUBECHIES_CENTERED,
	WL_HAAR,
	WL_HAAR_CENTERED,
	WL_BSPLINE,
	WL_BSPLINE_CENTERED

};


/****************
 ** DWT traits **
 ****************/
# include "DWTTraits.hpp"



/**
 * @brief 2D Discrete wavelet transform for Matrix template (from GSL)
 */
template <class T>
class DWT {

	typedef typename DWTTraits<T>::Type Type;

public:

	DWT (const size_t&, const wlfamily&, const size_t&);
	
	virtual ~DWT();
    
    Matrix<T> Trafo        (const Matrix<T>& m) const;
	Matrix<T> Adjoint      (const Matrix<T>& m) const;
    Matrix<T> operator*    (const Matrix<T>& m) const;
    Matrix<T> operator->* (const Matrix<T>& m) const;
	
};


# endif // __DWT_HPP__
