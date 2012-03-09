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
	DFT         (const Matrix<int> size, const Matrix<S> mask = Matrix<S>(1), 
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

	Matrix<size_t> m_size;
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

#endif
