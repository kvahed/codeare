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
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */
class DFT {
	
public:
	
	DFT         (const Matrix<int> size, const Matrix<double> mask = Matrix<double>(1), const Matrix<double> pc = Matrix<double>(1));
	
	virtual 
	~DFT        ();
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	Trafo       (Matrix<T>& m) const ;
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	Adjoint     (Matrix<T>& m) const;
	
	
private:
	
	/**
	 * Static class
	 */
	DFT()  {};
	
	
	bool m_initialised;

	Matrix<size_t> m_size;
	Matrix<double> m_mask;
	Matrix<double> m_pc;

	fftwf_plan     m_fwdplanf;
	fftwf_plan     m_bwdplanf;
	fftw_plan      m_fwdplan;
	fftw_plan      m_bwdplan;

	size_t         m_N;

	bool m_have_mask;
	bool m_have_pc;

	fftwf_complex* m_in;
	fftwf_complex* m_out;
};

#endif
