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
 *        This operator acts of Cartesian k-space data
 */
template <class T>
class CSENSE : public FT<T> {
	
public:
	

	/**
	 * @brief          Default constructor
	 */
	CSENSE() : m_initialised (false) {};


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

		m_dft  = new DFT<T> (ftdims);
		
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

		if (m_initialised)
			delete m_dft;
	
	}


	/**
	 * @brief          Forward transform
	 *
	 * @param  m       To transform
	 * @return         Transform
	 */
	Matrix< std::complex<T> >
	Adjoint       (const Matrix< std::complex<T> >& m) const {
		
		Matrix< std::complex<T> > res;
		Matrix< std::complex<T> > s (m_nc, m_af);
		Matrix< std::complex<T> > rp (m_nc);
		Matrix< std::complex<T> > tmp = m;

		// FT individual channels
		for (size_t i = 0; i < m_nc; i++)
			if (m_ndim == 2)
				if (m_af % 2)
					Slice  (tmp, i, *m_dft ->* fftshift(Slice  (tmp, i)));
				else
					Slice  (tmp, i, *m_dft ->*          Slice  (tmp, i));
			else
				if (m_af % 2)
					Volume (tmp, i, *m_dft ->* fftshift(Volume (tmp, i)));
				else
					Volume (tmp, i, *m_dft ->*          Volume (tmp, i));

		// Antialias
		for (size_t x = 0; x < m_dims[0]; x++)
			for (size_t y = 0; y < m_dims[1]/m_af; y++) 
				/*for (size_t z = 0; z < m_dims[2]; z++)*/ {

				for (size_t c = 0; c < m_nc; c++) {
					rp [c] = m (y, x, c);
					for (size_t i = 0; i < m_af; i++) 
						s (c, i) = m_sens (x, y + m_dims[1]/m_af * i, /*z,*/ c);
				}
				
				//s = inv(s'*s)*s';
				s = s.prodt(s);
				s = inv (s);                
				s = s.prodt(s);
				
				// rp=si*imfold(:,y,x);
				//rp = s.prodt(rp);
				/*	
						for (size_t i = 0; i < m_af; i++)
						res (y + m_dims[1] * i, x, z) = rp [i]; 
					*/
					
				}

		return tmp;
		
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
		
		assert (false);


	}


private:

	DFT<T>*        m_dft;

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
