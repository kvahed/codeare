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

#ifndef __CGRAPPA_HPP__
#define __CGRAPPA_HPP__

#include "FT.hpp"

/*
 * @brief GRAPPA operator<br/>
 *        Griswold et al. MRM 2002, vol. 47 (6) pp. 1202-1210
 */
template <class T>
class CGRAPPA : public FT<T> {



public:


	/**
	 * @brief    Default constructor
	 */
	CGRAPPA() {};


	/**
	 * @brief    Clean up and destroy
	 */
	virtual
	~CGRAPPA () {};


	/**
	 * @brief    Adjoint transform
	 */
	Matrix< std::complex<T> >
	Adjoint (Matrix< std::complex<T> >) const {
		Matrix<T> res;

		return res;
	};


	/**
	 * @brief    Adjoint transform
	 */
	Matrix< std::complex<T> >
	Trafo (Matrix< std::complex<T> >) const {
		Matrix<T> res;

		return res;
	};

private:


	Matrix< std::complex<T> > m_patch; /**< @brief Correction patch */
	Matrix< std::complex<T> > m_acs;   /**< @brief ACS lines        */

	DFT<T>**                  m_dft;   /**< @brief DFT operator     */

};

#endif /* __CGRAPPA_HPP__ */
