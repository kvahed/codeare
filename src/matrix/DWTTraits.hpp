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


/*****************
 ** gsl headers **
 *****************/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet2d.h>



/************************
 ** main struct (base) **
 ************************/
template <class T>
struct DWTTraits { };



/**
 * @brief C++ friendly interface to DWT (double precision)
 */
template<>
struct DWTTraits<double> {

	/**********************
	 ** type definitions **
	 **********************/
	typedef double mType;
	typedef double Type;
	
	
	// Allocate memory for type "Type"
	static inline
	Type *
	Malloc (size_t n) {
	  return (Type *) malloc ( n * sizeof (Type) );
	}
	
	
	// copy matrix data to real / imaginary array
	static inline
	void
	Prepare (Matrix<mType> & m, Type* re, Type* im) {
	  
   	memcpy (re, &m[0], m.Size() * sizeof(Type));
	  
	}
	
	
	// copy DWT data back to matrix
	static inline
	void
	Finalize (Matrix<mType> & m, Type* re, Type* im) {
	
    memcpy (&m[0], re, m.Size() * sizeof(Type));
	
	}
	
	
  // DWT non-standard, inverse
  static inline
  int
  NSTransformInverse (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
  
    if (!im)
      return gsl_wavelet2d_nstransform_inverse (w, data, tda, size1, size2, work);
    else
      return GSL_SUCCESS;
    
  }
	

  // DWT non-standard, forward
  static inline
  int
  NSTransformForward (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
		
		if (!im)
      return gsl_wavelet2d_nstransform_forward (w, data, tda, size1, size2, work);
    else
      return GSL_SUCCESS;

  }


};



/**
 * @brief C++ friendly interface to DWT (single precision)
 */
template<>
struct DWTTraits<float> {

	/**********************
	 ** type definitions **
	 **********************/
	typedef float mType;
	typedef double Type;
	
	
	// Allocate memory for type "Type"
	static inline
	Type *
	Malloc (size_t n) {
	  return (Type *) malloc ( n * sizeof (Type) );
	}
	
	
	// copy matrix data to real / imaginary array
	static inline
	void
	Prepare (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	  
   	for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++)
				re [i*sl+j] = m [j*sl+i];
	  
	}
	
	
	// copy DWT data back to matrix
	static inline
	void
	Finalize (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	
 		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) 
				m.At(j*sl+i) = re[i*sl+j];
	
	}
	
	
  // DWT non-standard, inverse
  static inline
  int
  NSTransformInverse (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
  
    if (!im)
      return gsl_wavelet2d_nstransform_inverse (w, data, tda, size1, size2, work);
    else
      return GSL_SUCCESS;
    
  }
	

  // DWT non-standard, forward
  static inline
  int
  NSTransformForward (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
		
		if (!im)
      return gsl_wavelet2d_nstransform_forward (w, data, tda, size1, size2, work);
    else
      return GSL_SUCCESS;

  }


};



/**
 * @brief C++ friendly interface to complex DWT (double precision)
 */
template<>
struct DWTTraits<cxdb> {


	/**********************
	 ** type definitions **
	 **********************/
	typedef cxdb mType;
	typedef double Type;
	
	
	// Allocate memory for type "Type"
	static inline
	Type *
	Malloc (size_t n) {
	  return (Type *) malloc ( n * sizeof (Type) );
	}


	// copy matrix data to real / imaginary array
	static inline
	void
	Prepare (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	  
 		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) {
				re [i*sl+j] = m.At(j*sl+i).real();
				im [i*sl+j] = m.At(j*sl+i).imag();
			}   	
   	
	}
	
	
	// copy DWT data back to matrix
	static inline
	void
	Finalize (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	
		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) 
				m.At(j*sl+i) = mType (re[i*sl+j],im [i*sl+j]);
 		
	}

	
  // complex DWT non-standard, inverse
  static inline
  int
  NSTransformInverse (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {

		return gsl_wavelet2d_nstransform_inverse (w, data, tda, size1, size2, work);
		
  }
	

  // complex DWT non-standard, forward
  static inline
  int
  NSTransformForward (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
		
    return gsl_wavelet2d_nstransform_forward (w, data, tda, size1, size2, work);
    
  }


};



/**
 * @brief C++ friendly interface to complex DWT (single precision)
 */
template<>
struct DWTTraits<cxfl> {


	/**********************
	 ** type definitions **
	 **********************/
	typedef cxfl mType;
	typedef double Type;
	
	
	// Allocate memory for type "Type"
	static inline
	Type *
	Malloc (size_t n) {
	  return (Type *) malloc ( n * sizeof (Type) );
	}


	// copy matrix data to real / imaginary array
	static inline
	void
	Prepare (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	  
 		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) {
				re [i*sl+j] = m.At(j*sl+i).real();
				im [i*sl+j] = m.At(j*sl+i).imag();
			}   	
   	
	}
	
	
	// copy DWT data back to matrix
	static inline
	void
	Finalize (Matrix<mType> & m, Type* re, Type* im, size_t sl) {
	
		for (size_t j = 0; j < sl; j++)
			for (size_t i = 0; i < sl; i++) 
				m.At(j*sl+i) = mType (re[i*sl+j],im [i*sl+j]);
 		
	}

	
  // complex DWT non-standard, inverse
  static inline
  int
  NSTransformInverse (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {

		return gsl_wavelet2d_nstransform_inverse (w, data, tda, size1, size2, work);
		
  }
	

  // complex DWT non-standard, forward
  static inline
  int
  NSTransformForward (const gsl_wavelet * w, Type * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work, bool im = false)
  {
		
    return gsl_wavelet2d_nstransform_forward (w, data, tda, size1, size2, work);
    
  }


};



# endif // __DWT_TRAITS_HPP__
