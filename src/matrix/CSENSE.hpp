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

#ifndef __CCSENSE_HPP__
#define __CCSENSE_HPP__

#include "Algos.hpp"
#include "CX.hpp"
#include "DFT.hpp"
#include "Access.hpp"

/**
 * @brief SENSE: Sensitivity Encoding for Fast MRI<br/>
 *        MRM (1999): vol. 42 (5) pp. 952-962<br/>
 *        This multi-threaded operator acts on Cartesian 
 *        k-space data 
 *
 * 
 */
template <class T>
class CSENSE : public FT<T> {
	
public:
	

	/**
	 * @brief          Default constructor
	 */
	CSENSE() {

		// Some initing
		m_initialised = false;
		m_af          = 1;
		m_dft         = 0;
		m_ndim        = 1;
		m_nc          = 1;

	}


	/**
	 * @brief          Construct CSENSE plans for forward and backward transform with credentials<br/>
	 *                 As of now we expect acceleration only in one direction
	 * 
	 * @param  sens    Sensitivity maps if imsize
	 * @param  af      Acceleration factor vector 2/3 elements for 2D/3D 
	 * @param  mask    K-Space mask
	 * @param  pc      Phase correction applied before forward or after adjoint transforms (default: empty)
	 * @param  b0      Off-resonance maps if available (default empty)
	 */
	CSENSE             (const Matrix< std::complex<T> >& sens, const unsigned short& af,
			            const Matrix<T>& mask = Matrix<T>(1), 
						const Matrix< std::complex<T> >& pc = Matrix< std::complex<T> >(1),
						const Matrix<T>& b0 = Matrix<T>(1)) : m_initialised (false) {

		// Some initing
		m_initialised = false;
		m_af          = 1;
		m_dft         = 0;
		m_ndim        = 1;
		m_nc          = 1;

		// Sensitivity maps dictate FT size
		m_dims = size(sens);

		// We expect sensitivities O (Ch,X,Y[,Z])
		m_ndim = numel(m_dims)-1;
		m_nc   = m_dims[m_ndim];
		
		// Handle only 2D / 3D
		assert (m_ndim == 2 || m_ndim == 3);
		
		// Set up FT
		Matrix<size_t> ftdims (m_ndim,1);
		
		for (size_t i = 0; i < m_ndim; i++)
			ftdims[i] = m_dims[i];
		
		ftdims[1] /= af;
		
		int np;

#pragma omp parallel default (shared)
		{
			np = omp_get_num_threads ();
		}	

		m_dft = new DFT<T>* [np];
		
		for (size_t i = 0; i < np; i++)
			m_dft[i]  = new DFT<T> (ftdims, mask, pc, b0);
		
		// Privates
		m_sens = sens;
		m_af   = af;
		
		// We're good
		m_initialised = true;
		
	}


	/**
	 * @brief          Clean up and destruct
	 */
	~CSENSE            () {

		int np;
		


#pragma omp parallel default (shared)
		{
			np = omp_get_num_threads ();
		}	
		


		if (m_initialised)
			for (int i = 0; i < np; i++)
				delete m_dft[i];
	
	}


	/**
	 * @brief          Forward transform
	 *
	 * @param  m       To transform
	 * @return         Transform
	 */
	Matrix< std::complex<T> >
	Adjoint       (const Matrix< std::complex<T> >& m) const {
		
		Matrix< std::complex<T> > res (m_dims[0], m_dims[1]/*, m_dims[2]*/);
		Matrix< std::complex<T> > tmp = m;
		
#pragma omp parallel
		{
			
			int tid = omp_get_thread_num ();
			
			Matrix< std::complex<T> > s  (m_nc, m_af);
			Matrix< std::complex<T> > si (m_af, m_af);
			Matrix< std::complex<T> > ra (m_nc,1);
			Matrix< std::complex<T> > rp (m_af,1);
			
#pragma omp for 
			
			// FT individual channels
			for (size_t i = 0; i < m_nc; i++)
				if (m_ndim == 2)
					Slice  (tmp, i, *(m_dft[tid]) ->* fftshift(Slice  (tmp, i)));
				else
					Volume (tmp, i, *(m_dft[tid]) ->* fftshift(Volume (tmp, i)));
			
			
#pragma omp for schedule (guided)
			
			// Antialias
			for (size_t x = 0; x < m_dims[0]; x++)
				for (size_t y = 0; y < m_dims[1]/m_af; y++) 
					/*for (size_t z = 0; z < m_dims[2]; z++)*/ {
					
					for (size_t c = 0; c < m_nc; c++) {
						
						ra [c] = tmp (x, y, c);
						
						for (size_t i = 0; i < m_af; i++) 
							s (c, i) = m_sens (x, y + m_dims[1]/m_af * i, /*z,*/ c);
						
					}
					
					//s = inv(s'*s)*s';
					si = gemm (s,   s, 'C', 'N');
					si = inv  (si);                
					si = gemm (si,  s, 'N', 'C');
					
					// rp=si*imfold(:,y,x);
					rp = gemm (si, ra, 'N', 'N');
					
					for (size_t i = 0; i < m_af; i++)
						res (x, y + m_dims[1]/m_af * i /*, z*/) = rp [i]; 
					
				}
			
		}		
		
		return res;
		
	}
	
	
	/**
	 * @brief          Backward transform (SENSE backward trafo? I don't know!)<br/> Bloedsinn!
	 *                 Why would anyone want to go back to sensitivity weighted undersampled k-space data?
	 *
	 * @param  m       To transform
	 * @return         Bummer! (This is not nice, right?)
	 */
	Matrix< std::complex<T> > 
	Trafo             (const Matrix< std::complex<T> >& m) const {
		
		Matrix < complex<T> > res;

		return res;


	}


private:

	DFT<T>**       m_dft;

	Matrix<T>      m_b0;

	Matrix< std::complex<T> > m_sens;
	Matrix< std::complex<T> > m_pc;
	
	unsigned short m_af;
	
	size_t         m_ndim;
	Matrix<size_t> m_dims; /**< Operator dimensionality Valid: [2,3]*/
	size_t         m_nc;

	bool           m_initialised;

};


#endif
