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

#include <fftw3.h>
#include "Matrix.hpp"
#include "Algos.hpp"
#include "IO.hpp"
#include "FT.hpp"


template <class T>
struct FTTraits { };


/**
 * @brief C++ friendly interface to complex FFTW (single precision)
 */
template<>
struct FTTraits<float> {
	
	typedef fftwf_plan    Plan; /**< @brief fftw plan (float precision) */
	typedef fftwf_complex Type; /**< @brief fftw complex data type (float precision) */
	

	/**
	 * @brief         DFT plan
	 *
	 * @param  rank   FT dimesionality
	 * @param  n      Size lengths of individual dimensions
	 * @param  in     Input memory
	 * @param  out    Output memory
	 * @param  sign   FT direction
	 * @param  flags  FFTW flags
	 *
	 * @return        Plan
	 */
	static inline Plan 
	DFTPlan (int rank, const int *n, Type *in, Type *out, int sign, unsigned flags) {                    
	    return fftwf_plan_dft (rank, n, in, out, sign, flags);
	}
	

	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline Type*
	Malloc (size_t n) {
		return fftwf_alloc_complex (n);
	}
	

	/**
	 * @brief        Free allocated memory
	 *
	 * @param  p     Memory address
	 */
	static inline void 
	Free (void* p) {
		fftwf_free (p);
	}
	

	/**
	 * @brief        Clean up tools
	 */
	static inline void 
	CleanUp () {
		fftwf_cleanup ();
	}
	

	/**
	 * @brief        Destroy plan
	 *
	 * @param        Plan to be destroyed
	 */
	static inline void 
	Destroy (Plan p) { 
		fftwf_destroy_plan (p); 
	}
	

	/**
	 * @brief       Execute plan
	 *
	 * @param       Plan to be executed
	 */
	static inline void 
	Execute (Plan p) { 
		fftwf_execute (p); 
	}        
	
};


/**
 * @brief C++ friendly interface to complex FFTW (double precision)
 */
template<>
struct FTTraits<double> {
	
	typedef fftw_plan    Plan;  /**< @brief fftw plan (double precision) */
	typedef fftw_complex Type;  /**< @brief fftw complex data type (double precision) */
	

	/**
	 * @brief         DFT plan
	 *
	 * @param  rank   FT dimesionality
	 * @param  n      Size lengths of individual dimensions
	 * @param  in     Input memory
	 * @param  out    Output memory
	 * @param  sign   FT direction
	 * @param  flags  FFTW flags
	 *
	 * @return        Plan
	 */
	static inline Plan 
	DFTPlan (int rank, const int *n, Type *in, Type *out, int sign, unsigned flags) {                    
	    return fftw_plan_dft (rank, n, in, out, sign, flags);
	}
	

	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline Type*
	Malloc (size_t n) {
		return fftw_alloc_complex (n);
	}
	

	/**
	 * @brief        Free allocated memory
	 *
	 * @param  p     Memory address
	 */
	static inline void
	Free (void* p) {
		fftw_free(p);
	}


	/**
	 * @brief        Clean up tools
	 */
	static inline void 
	CleanUp () {
		fftw_cleanup ();
	}
	

	/**
	 * @brief        Destroy plan
	 *
	 * @param   Plan Plan to be destroyed
	 */
	static inline void 
	Destroy (Plan p) { 
		fftw_destroy_plan (p); 
	}


	/**
	 * @brief       Execute plan
	 *
	 * @param  Plan Plan to be executed
	 */
	static inline void
	Execute (Plan p) { 
		fftw_execute (p); 
	}        

};


/**
 * @brief         FFT shift
 * 
 * @param   m     TO be shifted
 * @return        Shifted
 */
template<class T> inline static Matrix<T>
fftshift (const Matrix<T>& m) {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
	Matrix<T> res  = m;
	
	for (size_t s = 0; s < m.Dim(2); s++)
		for (size_t l = 0; l < m.Dim(1); l++)
			for (size_t c = 0; c < m.Dim(0); c++)
				res (c,l,s) *= (T) pow ((T)-1.0, (T)(s+l+c));
	
	return res;
	
}
	

/**
 * @brief         Hann window
 * 
 * @param   size  Side lengths
 * @param   t     Scaling factor
 * @return        Window
 */
template <class T> inline static Matrix< std::complex<T> >
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
				 const Matrix< std::complex<T> >& pc = Matrix< std::complex<T> >(1),
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
				 const Matrix< std::complex<T> >& pc = Matrix< std::complex<T> >(1),
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
	

	/**
	 * @brief        Clean up RAM, destroy plans
	 */
	virtual 
	~DFT        () {

		FTTraits<T>::Destroy (m_fwplan);
		FTTraits<T>::Destroy (m_bwplan);
		
		FTTraits<T>::CleanUp();
		
		FTTraits<T>::Free (m_in); 
		FTTraits<T>::Free (m_out);
		
	}
	
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	Trafo       (const Matrix< std::complex<T> >& m) const {
		
		Matrix< std::complex<T> > res = fftshift (m);
		
		memcpy (m_in, &res[0], m_cs);
		
		if (m_have_pc)
			res *= m_pc;
		
		FTTraits<T>::Execute (m_fwplan);
		
		memcpy (&res[0], m_out, m_cs);
		
		if (m_have_mask)
			res *= m_mask;

		res  = fftshift (res) / m_sn;
		
		return res;
		
	}
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	Adjoint     (const Matrix< std::complex<T> >& m) const {
		
		Matrix< std::complex<T> > res = fftshift (m);
		
		if (m_have_mask)
			res *= m_mask;
		
		memcpy (m_in,  &res[0], m_cs);
		
		FTTraits<T>::Execute (m_bwplan);
		
		memcpy (&res[0], m_out, m_cs);
		
		if (m_have_pc)
			res *= m_cpc;
		
		res  = fftshift (res) / m_sn;

		return res;
			
	}
	
	
private:
	

	/**
	 * @brief       Allocate RAM and plans
	 *
	 * @param  rank FT rank
	 * @param  n    Side lengths
	 */
	void 
	Allocate (const int rank, const int* n) {
		
		m_in     = FTTraits<T>::Malloc  (m_N);
		m_out    = FTTraits<T>::Malloc  (m_N);
		m_fwplan = FTTraits<T>::DFTPlan (rank, n, m_in, m_out, FFTW_FORWARD,  FFTW_MEASURE);
		m_bwplan = FTTraits<T>::DFTPlan (rank, n, m_in, m_out, FFTW_BACKWARD, FFTW_MEASURE);

		m_cs     = m_N * sizeof(Type);
		m_sn     = sqrt (m_N);

	}


	bool      m_initialised;  /**< @brief Memory allocated / Plans, well, planned! :)*/
	
	Matrix<T> m_mask;         /**< @brief K-space mask (applied before inverse and after forward transforms) (double precision)*/
	
	Matrix< std::complex<T> > 
	          m_pc;           /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	Matrix< std::complex<T> > 
	          m_cpc;          /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	
	Plan      m_fwplan;       /**< @brief Forward plan (double precision)*/
	Plan      m_bwplan;       /**< @brief Backward plan (double precision)*/

	size_t    m_N;            /**< @brief # Nodes */
	size_t    m_cs;

	T         m_sn;

	bool      m_have_mask;    /**< @brief Apply mask?*/
	bool      m_have_pc;      /**< @brief Apply phase correction?*/
	bool      m_zpad;         /**< @brief Zero padding? (!!!NOT OPERATIONAL YET!!!)*/
	
	Type*     m_in;           /**< @brief Aligned fftw input*/
	Type*     m_out;          /**< @brief Aligned fftw output*/

	
};



#endif




