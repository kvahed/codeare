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


#ifndef __DWT_TRAITS_HPP__
#define __DWT_TRAITS_HPP__


#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet2d.h>

#include <string.h>

template <class T>
struct DWTTraits;


/**
 * @brief C++ friendly interface to complex DWT (double precision)
 */
template<>
struct DWTTraits<cxfl> {

	typedef cxfl mType;
	typedef double Type;

	// copy matrix data to real / imaginary array
	static inline void
	prepare (const Matrix<mType>& m, VECTOR_TYPE(Type)& re, VECTOR_TYPE(Type)& im, const size_t& sl) {
 		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) {
				re [i*sl+j] = m.At(j*sl+i).real();
				im [i*sl+j] = m.At(j*sl+i).imag();
			}
	}


	// copy DWT data back to matrix
	static inline void
	finalize (Matrix<mType>& m, const VECTOR_TYPE(Type)& re, const VECTOR_TYPE(Type)& im, const size_t& sl) {
		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++)
				m[j*sl+i] = mType (re[i*sl+j], im[i*sl+j]);
	}


	// complex DWT non-standard, inverse
	static inline int
	nstransform_inverse (const gsl_wavelet* w, VECTOR_TYPE(Type)& data, const size_t& tda,
			const size_t& size1, const size_t& size2, gsl_wavelet_workspace* work, const bool& im = false) {
		return gsl_wavelet2d_nstransform_inverse (w, &data[0], tda, size1, size2, work);
	}


	// complex DWT non-standard, forward
	static inline int
	nstransform_forward (const gsl_wavelet* w, VECTOR_TYPE(Type)& data, const size_t& tda,
			const size_t& size1, const size_t& size2, gsl_wavelet_workspace* work, const bool& im = false) {
		return gsl_wavelet2d_nstransform_forward (w, &data[0], tda, size1, size2, work);
	}


};


/**
 * @brief C++ friendly interface to complex DWT (double precision)
 */
template<>
struct DWTTraits<cxdb> {


	typedef cxdb mType;
	typedef double Type;
	
	// copy matrix data to real / imaginary array
	static inline void
	prepare (const Matrix<mType>& m, VECTOR_TYPE(Type)& re, VECTOR_TYPE(Type)& im, const size_t& sl) {
 		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) {
				re [i*sl+j] = m.At(j*sl+i).real();
				im [i*sl+j] = m.At(j*sl+i).imag();
			}   	
	}
	
	
	// copy DWT data back to matrix
	static inline void
	finalize (Matrix<mType>& m, const VECTOR_TYPE(Type)& re, const VECTOR_TYPE(Type)& im, const size_t& sl) {
		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++)
				m[j*sl+i] = mType (re[i*sl+j], im[i*sl+j]);
	}

	
	// complex DWT non-standard, inverse
	static inline int
	nstransform_inverse (const gsl_wavelet* w, VECTOR_TYPE(Type)& data, const size_t& tda,
			const size_t& size1, const size_t& size2, gsl_wavelet_workspace* work, const bool& im = false) {
		return gsl_wavelet2d_nstransform_inverse (w, &data[0], tda, size1, size2, work);
	}


	// complex DWT non-standard, forward
	static inline int
	nstransform_forward (const gsl_wavelet* w, VECTOR_TYPE(Type)& data, const size_t& tda,
			const size_t& size1, const size_t& size2, gsl_wavelet_workspace* work, const bool& im = false) {
		return gsl_wavelet2d_nstransform_forward (w, &data[0], tda, size1, size2, work);
	}
	

};


# endif // __DWT_TRAITS_HPP__
