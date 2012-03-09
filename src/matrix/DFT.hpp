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

#ifndef __DFT_HPP__
#define __DFT_HPP__

#include <fftw3.h>
#include "Matrix.hpp"
#include "Algos.hpp"

/**
 * @brief Matrix templated 1-3D Discrete Cartesian Fourier transform
 */
template <class T>
class DFT {
	
public:
	

	/**
	 * @brief        Construct FFTW plans for forward and backward FT with credentials
	 * 
	 * @param  size  Matrix of side length of the FT range
	 * @param  mask  K-Space mask (if left empty no mask is applied)
	 * @param  pc    Phase correction (or target phase)
	 */
	template<class S>
	DFT         (const Matrix<size_t>& size, const Matrix<S> mask = Matrix<S>(1), 
				 const Matrix<T> pc = Matrix<T>(1));
	

	/**
	 * @brief        Construct FFTW plans for forward and backward FT with credentials for FT with identical side lengths
	 * 
	 * @param  rank  Rank (i.e. # FT directions)
	 * @param  sl    Side length of the slice, volume ...
	 * @param  mask  K-Space mask (if left empty no mask is applied)
	 * @param  pc    Phase correction (or target phase)
	 */
	template<class S>
	DFT         (const size_t rank, const size_t sl, const Matrix<S> mask = Matrix<S>(1), 
				 const Matrix<T> pc = Matrix<T>(1));
	

	/**
	 * @brief        Construct FFTW plans for forward and backward FT with credentials
	 * 
	 * @param  size  Matrix of side length of the FT range
	 * @param  mask  K-Space mask (if left empty no mask is applied)
	 * @param  pc    Phase correction (or target phase)
	 */
	DFT         (const Matrix<size_t>& size);
	

	/**
	 * @brief        Clean up RAM, destroy plans
	 */
	virtual 
	~DFT        ();
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Trafo       (const Matrix<T>& m) const ;
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Adjoint     (const Matrix<T>& m) const;
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	operator*   (const Matrix<T>& m) const {

		return Trafo(m);

	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	operator->* (const Matrix<T>& m) const {
		
		return Adjoint (m);

	}
	

	template <class S> 
	void Mask (const Matrix<S>& m);

private:
	
	bool m_initialised;

	Matrix<double> m_mask;
	Matrix<T>      m_pc;
	Matrix<T>      m_cpc;
	
	fftwf_plan     m_fwdplanf;
	fftwf_plan     m_bwdplanf;
	fftw_plan      m_fwdplan;
	fftw_plan      m_bwdplan;
	
	size_t         m_N;
	
	bool           m_have_mask;
	bool           m_have_pc;
	
	void*          m_in;
	void*          m_out;

};

template <class T> inline static Matrix<cxfl> 
fftshift        (const Matrix<T>& m) {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
	Matrix<T> res  = m;
	
	for (size_t s = 0; s < m.Dim(2); s++)
		for (size_t l = 0; l < m.Dim(1); l++)
			for (size_t c = 0; c < m.Dim(0); c++)
				res.At (c,l,s) *= (float) pow ((float)-1.0, (float)(s+l+c));
		
	return res;
	
}


template <class T> inline static Matrix<cxfl> 
ifftshift        (const Matrix<T>& m) {
	
	return fftshift (m);	
	
}


inline Matrix<double> 
hannwindow (const Matrix<size_t>& size) {

	size_t dim = size.Dim(0);
	
	assert (dim > 1 && dim < 4);
	
	Matrix<double> res;
	
	if      (dim == 1) res = Matrix<double> (size[0], 1);
	else if (dim == 2) res = Matrix<double> (size[0], size[1]);
	else               res = Matrix<double> (size[0], size[1], size[2]);
	
	float          h, d;
	float          m[3];
	
	if (Is1D(res)) {
		
		m[0] = 0.5 * size[0];
		m[1] = 0.0;
		m[2] = 0.0;
		
	} else if (Is2D(res)) {
		
		m[0] = 0.5 * size[0];
		m[1] = 0.5 * size[1];
		m[2] = 0.0;
		
	} else {
		
		m[0] = 0.5 * size[0];
		m[1] = 0.5 * size[1];
		m[2] = 0.5 * size[2];
		
	}
	
	res = Squeeze(res);
	
	for (size_t s = 0; s < res.Dim(2); s++)
		for (size_t r = 0; r < res.Dim(1); r++)
			for (size_t c = 0; c < res.Dim(0); c++) {
				d = pow( (float)pow(((float)c-m[0])/m[0],2) + pow(((float)r-m[1])/m[1],2) + pow(((float)s-m[2])/m[2],2) , (float)0.5);
				h = (d < 1) ? (0.5 + 0.5 * cos (PI * d)) : 0.0;
				res(c,r,s) = h;
			}
	
	return res;
	
}

#endif




