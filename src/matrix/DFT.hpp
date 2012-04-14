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
#include "IO.hpp"
#include "FT.hpp"


/**
 * @brief Matrix templated 1-3D Discrete Cartesian Fourier transform
 */
template <class T>
class DFT : public FT<T> {
	
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
	 * @brief        Construct FFTW plans for forward and backward FT
	 * 
	 * @param  size  Matrix of side length of the FT range
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
	virtual Matrix<T> 
	Trafo       (const Matrix<T>& m) const ;
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> 
	Adjoint     (const Matrix<T>& m) const;
	
	
private:
	
	bool           m_initialised;  /**< @brief Memory allocated / Plans, well, planned! :)*/

	Matrix<double> m_mask;         /**< @brief K-space mask (applied before inverse and after forward transforms) (double precision)*/
	Matrix<float>  m_maskf;        /**< @brief K-space mask (applied before inverse and after forward transforms) (float precision)*/

	Matrix<T>      m_pc;           /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	Matrix<T>      m_cpc;          /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	
	fftwf_plan     m_fwdplanf;     /**< @brief Forward plan (double precision)*/
	fftwf_plan     m_bwdplanf;     /**< @brief Backward plan (double precision)*/
	fftw_plan      m_fwdplan;      /**< @brief Forward plan (single precision)*/
	fftw_plan      m_bwdplan;      /**< @brief Forward plan (single precision)*/
	
	size_t         m_N;            /**< @brief # Nodes */
	
	bool           m_have_mask;    /**< @brief Apply mask?*/
	bool           m_have_pc;      /**< @brief Apply phase correction?*/
	bool           m_zpad;         /**< @brief Zero padding? (!!!NOT OPERATIONAL YET!!!)*/
	
	void*          m_in;           /**< @brief Aligned fftw input*/
	void*          m_out;          /**< @brief Aligned fftw output*/

};


/**
 * @brief         FFT shift
 * 
 * @param   m     TO be shifted
 * @return        Shifted
 */
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


/**
 * @brief         FFT shift (Inverse: i.e. Forward - works only for matrices with even side lengths)
 * 
 * @param   m     TO be shifted
 * @return        Shifted
 */
template <class T> inline static Matrix<cxfl> 
ifftshift        (const Matrix<T>& m) {
	
	return fftshift (m);	
	
}


/**
 * @brief         Hann window
 * 
 * @param   size  Side lengths
 * @param   t     Scaling factor
 * @return        Window
 */
template <class T> inline static Matrix<T> 
hannwindow (const Matrix<size_t>& size, const T& t) {

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
	
	res = squeeze(res);
	
	for (size_t s = 0; s < res.Dim(2); s++)
		for (size_t r = 0; r < res.Dim(1); r++)
			for (size_t c = 0; c < res.Dim(0); c++) {
				d = pow( (float)pow(((float)c-m[0])/m[0],2) + pow(((float)r-m[1])/m[1],2) + pow(((float)s-m[2])/m[2],2) , (float)0.5);
				h = (d < 1) ? (0.5 + 0.5 * cos (PI * d)) : 0.0;
				res(c,r,s) = t * h;
			}
	
	return res;
	
}

template<> template<>
DFT<cxfl>::DFT (const size_t rank, const size_t sl, const Matrix<float> mask, const Matrix<cxfl> pc) :
	m_have_mask (false),
	m_have_pc (false) {
	
	int n[rank];
	
	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_maskf = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);

	for (size_t i = 0; i < rank; i++)
		n[i]  = sl;

	m_N   = pow (sl, rank);

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}

template<> template<>
DFT<cxfl>::DFT (const Matrix<size_t>& size, const Matrix<float> mask, const Matrix<cxfl> pc) : m_N(1),
																						m_have_mask (false),
																						m_have_pc (false) {

	int rank = size.Size();
	int n[rank];

	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_maskf = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (int i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_ESTIMATE);

	m_initialised = true;

}


template<> 
DFT<cxfl>::DFT (const Matrix<size_t>& size) : m_N(1),
											 m_have_mask (false),
											 m_have_pc (false) {

	
	int rank = size.Size();
	int n[rank];

	for (int i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}


template<> template<>
DFT<cxdb>::DFT (const size_t rank, const size_t sl, const Matrix<double> mask, const Matrix<cxdb> pc) :
	m_have_mask (false),
	m_have_pc (false) {
	
	int n[rank];
	
	if (numel(mask) > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (numel(pc) > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (size_t i = 0; i < rank; i++)
		n[i]  = sl;

	m_N   = pow (sl, rank);

	m_in  = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);
	m_out = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);

	m_fwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}

template<> template<>
DFT<cxdb>::DFT (const Matrix<size_t>& size, const Matrix<double> mask, const Matrix<cxdb> pc) : m_N(1),
																						m_have_mask (false),
																						m_have_pc (false) {

	int rank = size.Size();
	int n[rank];

	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (int i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);
	m_out = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);

	m_fwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_BACKWARD, FFTW_ESTIMATE);

	m_initialised = true;

}


template<>
DFT<cxfl>::~DFT () {

	fftwf_destroy_plan (m_fwdplanf);
	fftwf_destroy_plan (m_bwdplanf);

	fftwf_cleanup ();

	fftwf_free(m_in); 
	fftwf_free(m_out);

}


template<>
DFT<cxdb>::~DFT () {

	fftw_destroy_plan (m_fwdplan);
	fftw_destroy_plan (m_bwdplan);

	fftw_cleanup ();

	fftw_free(m_in); 
	fftw_free(m_out);

}


template<> Matrix<cxfl> 
DFT<cxfl>::Trafo (const Matrix<cxfl>& m) const {
	
    Matrix<cxfl> res = fftshift(m);
	memcpy (m_in, &res[0], sizeof(fftwf_complex) * m_N);
	if (m_have_pc)
		res *= m_pc;

	fftwf_execute(m_fwdplanf);

	memcpy (&res[0], m_out, sizeof(fftwf_complex) * m_N);
	if (m_have_mask)
		res *= m_maskf;

    return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxfl>
DFT<cxfl>::Adjoint (const Matrix<cxfl>& m) const {

    Matrix<cxfl> res = fftshift(m);
	if (m_have_mask)
		res *= m_maskf;
	memcpy (m_in, &res[0], sizeof(fftwf_complex) * m_N);

	fftwf_execute(m_bwdplanf);

	memcpy (&res[0], m_out, sizeof(fftwf_complex) * m_N);
	if (m_have_pc)
		res *= m_cpc;

	return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxdb> 
DFT<cxdb>::Trafo (const Matrix<cxdb>& m) const {
	
    Matrix<cxdb> res = fftshift(m);
	memcpy (m_in, &res[0], sizeof(fftw_complex) * m_N);
	if (m_have_pc)
		res *= m_pc;

	fftw_execute(m_fwdplan);

	memcpy (&res[0], m_out, sizeof(fftw_complex) * m_N);
	if (m_have_mask)
		res *= m_mask;

    return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxdb>
DFT<cxdb>::Adjoint (const Matrix<cxdb>& m) const {

    Matrix<cxdb> res = fftshift(m);
	if (m_have_mask)
		res *= m_mask;
	memcpy (m_in, &res[0], sizeof(fftw_complex) * m_N);

	fftw_execute(m_bwdplan);

	memcpy (&res[0], m_out, sizeof(fftw_complex) * m_N);
	if (m_have_pc)
		res *= m_cpc;

	return fftshift(res/sqrt((float)m.Size()));
	
}

#endif




