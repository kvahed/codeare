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
#include "Creators.hpp"

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
	CSENSE() : m_dft(0) {}


	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	CSENSE        (const Params& params) :
		FT<T>::FT(params), m_dft(0) {

		p = params;

		p.Set ("initialised", false);

		std::string map_name = p.Get<std::string> ("smaps_name");
		std::string img_name = p.Get<std::string> ("fimgs_name");

		sens = Workspace::Instance()->Get<std::complex<T> >(map_name);
		Matrix<std::complex<T> >& imgs = Workspace::Instance()->Get<std::complex<T> >(img_name);

		const size_t nc = size (sens, ndims(sens)-1);
		assert (nc > 1);
		p.Set ("nc", nc);

		size_t af = size(sens, 1) / size(imgs, 1);
		p.Set("af", af);

		size_t ndim = ndims(sens)-1;
		assert (ndim == 2 || ndim == 3);
		p.Set ("ndim", ndim);

		// We expect sensitivities O (X,Y,Z,CH)
		dims = ones<size_t> (3,1);
		for (size_t i = 0; i < ndim; i++)
			dims[i] = size(sens,i);
		dims[1] /= af;

		p.Set ("dims", dims);

		if (!p.exists("treg"))
			p.Set("treg", (T)0.0);

		// Multi-threading will need multiple FFTW plans
		int np;

#pragma omp parallel default (shared)
		{
			np = omp_get_num_threads ();
		}

		Matrix<size_t> ftdims = resize(dims,ndim,1);

		m_dft = new DFT<T>* [np];

		for (size_t i = 0; i < np; i++)
			m_dft[i]  = new DFT<T> (dims/*, mask, pc, b0*/);

		// Great
		p.Set ("initialised", true);

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

		if (m_dft)
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
		
		bool compgfm     = p.Get<bool>("compgfm");
		size_t af        = p.Get<size_t>("af");
		size_t nc        = p.Get<size_t>("nc");
		size_t ndim      = p.Get<size_t>("ndim");
		T      treg      = p.Get<T>("treg");
		bool initialised = p.Get<bool>("initialised");

		Matrix< std::complex<T> > res (dims[0], dims[1]*af, dims[2], (compgfm) ? 2 : 1);
		Matrix< std::complex<T> > tmp = m;
		
#pragma omp parallel
		{
			
			int tid = omp_get_thread_num ();
			
			Matrix<std::complex<T> > s  (nc, af);
			Matrix<std::complex<T> > si (af, af);
			Matrix<std::complex<T> > ra (nc,  1);
			Matrix<std::complex<T> > rp (af,  1);
			Matrix<std::complex<T> > gf (af,  1);
			Matrix<std::complex<T> > reg = treg * eye<std::complex<T> >(af);
			
#pragma omp for 
			
			// FT individual channels
			for (size_t i = 0; i < nc; i++)
				if (ndim == 2)
					Slice  (tmp, i, *(m_dft[tid]) ->* Slice  (tmp, i));
				else
					Volume (tmp, i, *(m_dft[tid]) ->* Volume (tmp, i));
			
#pragma omp for schedule (guided)
			
			// Antialiasing
			for (size_t x = 0; x < dims[0]; x++)
				for (size_t y = 0; y < dims[1]; y++)
					for (size_t z = 0; z < dims[2]; z++) {
						
						for (size_t c = 0; c < nc; c++) {
							
							ra [c] = (dims[2]-1) ? tmp (x, y, z, c) : tmp (x, y, c);
							
							for (size_t i = 0; i < af; i++)
								s (c, i) = (dims[2]-1) ? sens (x, y + dims[1] * i, z, c) : sens (x, y + dims[1] * i, c);
							
						}
						
						si = gemm (s,   s, 'C', 'N');

						if (compgfm)
							gf = diag (si);

						si = inv  (si + reg);

						if (compgfm)
							gf = diag (si) * gf;

						si = gemm (si,  s, 'N', 'C');
						rp = gemm (si, ra, 'N', 'N');
						
						for (size_t i = 0; i < af; i++) {

							res (x, y + dims[1] * i, z, 0) =          rp [i];

							if (compgfm)
								res (x, y + dims[1] * i, z, 1) = sqrt(abs(gf [i]));

						}

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

	DFT<T>**              m_dft;
	Params                p;
	Matrix < complex<T> > sens;
	Matrix <size_t>       dims;


};


#endif
