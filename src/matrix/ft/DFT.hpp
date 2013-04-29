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

#ifndef __DFT_HPP__
#define __DFT_HPP__

#include "Matrix.hpp"
#include "Algos.hpp"
#include "IO.hpp"
#include "FT.hpp"
#include "FFTWTraits.hpp"


/**
 * @brief         FFT shift
 * 
 * @param   m     TO be shifted
 * @return        Shifted
 */
template<class T> inline Matrix<T>
fftshift (const Matrix<T>& m) {

	size_t cc = floor(((float)size(m,COL))/2);
	size_t cl = floor(((float)size(m,LIN))/2);

	Matrix<size_t> d = size(m);
	size_t ii,jj;

	Matrix<T> res = m;

	for (int i = 0; i < d[1]; i++) {
	    ii = (i + cl) % d[1];
	    for (int j = 0; j < d[0]; j++) {
	    	jj = (j + cc) % d[0];
	    	res(jj,ii) = m(j,i);
	    }
	}

	return res;

}
	

/**
 * @brief         Inverse FFT shift
 *
 * @param   m     TO be inversely shifted
 * @return        Inversely shifted
 */
template<class T> inline Matrix<T>
ifftshift (const Matrix<T>& m) {

	size_t cc = floor(((float)size(m,COL))/2);
	size_t cl = floor(((float)size(m,LIN))/2);

	Matrix<size_t> d = size(m);
	size_t ii,jj;

	Matrix<T> res = m;

	for (int i = 0; i < d[1]; i++) {
	    ii = (i + cl) % d[1];
	    for (int j = 0; j < d[0]; j++) {
	    	jj = (j + cc) % d[0];
	    	res(j,i) = m(jj,ii);
	    }
	}

	return res;

}


/**
 * @brief         Hann window
 * 
 * @param   size  Side lengths
 * @param   t     Scaling factor
 * @return        Window
 */
template <class T> inline Matrix< std::complex<T> >
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


/**
 * @brief Matrix templated 1-3D Discrete Cartesian Fourier transform
 */
template <class T>
class DFT : public FT<T> {


	typedef typename FTTraits<T>::Plan Plan;
	typedef typename FTTraits<T>::Type Type;
	typedef typename std::complex<T>   CT;


	
public:
	
	
	/**
	 * @brief        Construct FFTW plans for forward and backward FT with credentials
	 * 
	 * @param  sl    Matrix of side length of the FT range
	 * @param  mask  K-Space mask (if left empty no mask is applied)
	 * @param  pc    Phase correction (or target phase)
	 * @param  b0    Field distortion
	 */
	DFT         (const Matrix<size_t>& sl, const Matrix<T>& mask,
				 const Matrix<CT>& pc = Matrix<CT>(1),
				 const Matrix<T>& b0 = Matrix<T>(1)) :
		m_N(1), m_have_mask (false), m_have_pc (false) {

		size_t rank = numel(sl);

		int n [rank];
		
		if (numel(mask) > 1) {
			m_have_mask = true;
			m_mask      = mask;
		}
		
		if (numel(pc)   > 1) {
			m_have_pc   = true;
			m_pc        = pc;
			m_cpc       = conj(pc);
		}
		
		for (int i = 0; i < rank; i++) {
			n[i]  = (int) sl [rank-1-i];
			m_N  *= n[i];
		}
		
		Allocate (rank, n);
		
		m_initialised = true;

	}

	
	/**
	 * @brief        Construct FFTW plans for forward and backward FT with credentials for FT with identical side lengths
	 * 
	 * @param  rank  Rank (i.e. # FT directions)
	 * @param  sl    Side length of the slice, volume ...
	 * @param  mask  K-Space mask (if left empty no mask is applied)
	 * @param  pc    Phase correction (or target phase)
	 * @param  b0    Static field distortion
	 */
	DFT         (const size_t rank, const size_t sl, const Matrix<T>& mask = Matrix<T>(1),
				 const Matrix<CT>& pc = Matrix<CT>(1),
				 const Matrix<T>& b0 = Matrix<T>(1)) :
		m_have_mask (false), m_have_pc (false) {

		int n[rank];
		
		if (numel(mask) > 1) {
			m_have_mask = true;
			m_mask      = mask;
		}
		
		if (pc.Size() > 1) {
			m_have_pc   = true;
			m_pc   = pc;
			m_cpc  = conj(pc);
		}
		
		for (size_t i = 0; i < rank; i++) {
			n[i]  = sl;
			m_N  *= n[i];
		}
		
		Allocate (rank, n);
		
		m_initialised = true;
	
	}
	
	
	/**
	 * @brief        Construct FFTW plans for forward and backward FT
	 * 
	 * @param  sl    Matrix of side lengths of the FT range
	 */
	DFT         (const Matrix<size_t>& sl) : 
		m_N(1), m_have_mask (false), m_have_pc (false) {
		
		
		int rank = numel(sl);
		int n[rank];
		
		for (int i = 0; i < rank; i++) {
			n[i]  = (int) sl[rank-1-i];
			m_N  *= n[i];
		}

		Allocate (rank, n);
		
		m_initialised = true;
		
	}
	

	DFT        (const Params& params) :
		FT<T>::FT(params), m_cs(0), m_N(0), m_in(0), m_have_pc(false), m_zpad(false),
		m_initialised(false), m_have_mask(false) {

	}


    DFT () :
    	m_cs(0), m_N(0), m_in(0), m_have_pc(false), m_zpad(false),
    	m_initialised(false), m_have_mask(false){}
    
	/**
	 * @brief        Clean up RAM, destroy plans
	 */
	virtual 
	~DFT        () {

		FTTraits<T>::Destroy (m_fwplan);
		FTTraits<T>::Destroy (m_bwplan);
		
		FTTraits<T>::CleanUp();
		
		FTTraits<T>::Free (m_in);

	}
	
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline virtual Matrix<CT>
	Trafo       (const Matrix<CT>& m) const {
		
		Matrix<CT> res = ifftshift((m_have_pc) ? m * m_pc : m);

		FTTraits<T>::Execute (m_fwplan, (Type*)&res[0], (Type*)&res[0]);

		res = fftshift(res);

		if (m_have_mask)
			res *= m_mask;
		
		return res / m_sn;
		
	}
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline virtual Matrix<CT>
	Adjoint     (const Matrix<CT>& m) const {

		Matrix<CT> res = ifftshift((m_have_mask) ? m * m_mask : m);

		FTTraits<T>::Execute (m_bwplan, (Type*)&res[0], (Type*)&res[0]);

		res = fftshift(res);

		if (m_have_pc)
			res *= m_cpc;
		
		return res / m_sn;
			
	}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline virtual Matrix<CT>
	operator* (const Matrix<CT>& m) const {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline virtual Matrix<CT>
	operator->* (const Matrix<CT>& m) const {
		return Adjoint (m);
	}


private:
	

	/**
	 * @brief       Allocate RAM and plans
	 *
	 * @param  rank FT rank
	 * @param  n    Side lengths
	 */
	inline void
	Allocate (const int rank, const int* n) {
		
		m_in     = FTTraits<T>::Malloc  (m_N);

		m_fwplan = FTTraits<T>::DFTPlan (rank, n, m_in, m_in, FFTW_FORWARD,  FFTW_MEASURE);
		m_bwplan = FTTraits<T>::DFTPlan (rank, n, m_in, m_in, FFTW_BACKWARD, FFTW_MEASURE);

		m_cs     = m_N * sizeof(Type);
		m_sn     = sqrt (m_N);

	}


	bool       m_initialised;  /**< @brief Memory allocated / Plans, well, planned! :)*/
	
	Matrix<T>  m_mask;         /**< @brief K-space mask (applied before inverse and after forward transforms) (double precision)*/
	
	Matrix<CT> m_pc;           /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	Matrix<CT> m_cpc;          /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	
	Plan       m_fwplan;       /**< @brief Forward plan (double precision)*/
	Plan       m_bwplan;       /**< @brief Backward plan (double precision)*/

	size_t     m_N;            /**< @brief # Nodes */
	size_t     m_cs;

	T          m_sn;

	bool       m_have_mask;    /**< @brief Apply mask?*/
	bool       m_have_pc;      /**< @brief Apply phase correction?*/
	bool       m_zpad;         /**< @brief Zero padding? (!!!NOT OPERATIONAL YET!!!)*/
	
	Type*      m_in;           /**< @brief Aligned fftw input*/

	
};



#endif




