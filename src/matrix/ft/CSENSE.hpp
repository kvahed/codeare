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
#include "Lapack.hpp"

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
	
	typedef std::complex<T> CT;

	/**
	 * @brief          Default constructor
	 */
	CSENSE() : m_dft(0), nthreads(1), compgfm(0), aaf(1), nc(1), initialised(0), ndim(1) {}


	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	CSENSE        (const Params& params) :
		FT<T>::FT(params), m_dft(0), nthreads (1), treg(0.), compgfm(0) {

		Workspace& w = Workspace::Instance();

		// Maps & 1st set of images
		sens = w.Get<CT>(params.Get<std::string> ("smaps_name"));
		Matrix<CT>& imgs = w.Get<CT>(params.Get<std::string> ("fimgs_name"));

		// Channels
		nc = size (sens, ndims(sens)-1);
		assert (nc > 1);

		// 3D acceleration
		af = size(sens) / size(imgs);
		std::cout << "\n  acceleration vector: " << af;
		aaf = prod (af);


		// Can only handle 2D/3D SENSE
		ndim = ndims(sens)-1;
		assert (ndim == 2 || ndim == 3);

		// We expect sensitivities O (X,Y,Z,CH)
		dims = size(imgs);
		dims = resize(dims,1,size(dims,1)-1);
		std::cout << "  fft dims: " << dims ;

		TikhonovMat(params);

		AllocateDFTs(params);

		// Need 3-dimensional extents even if 2D
		if (numel(dims) == 2) {
			dims = resize (dims,1,3);
			dims[2] = 1;
		}

		// We're good
		initialised = true;

	}


	/**
	 * @brief          Clean up and destruct
	 */
	virtual
	~CSENSE            () {

		if (initialised)
			for (int i = 0; i < nthreads; i++)
				if (m_dft[i] != 0)
					delete m_dft[i];

	}


	/**
	 * @brief          Forward transform
	 *
	 * @param  m       To transform
	 * @return         Transform
	 */
	Matrix<CT>
	Adjoint       (const Matrix<CT>& m) const {

		Matrix<CT> res (dims[0]*af[0], dims[1]*af[1], (ndim == 3) ? dims[2]*af[2] : 1);
		Matrix<CT> tmp = m;

		omp_set_num_threads(nthreads);

#pragma omp parallel
		{
			
			int tid = omp_get_thread_num ();
			Matrix<CT> s  (nc,  aaf);
			Matrix<CT> si (aaf, aaf);
			Matrix<CT> ra (nc,    1);
			Matrix<CT> rp (aaf,   1);
			Matrix<CT> gf (aaf,   1);
			
			DFT<T>& ft = *(m_dft[tid]);

#pragma omp for
			
			// FT individual channels
			for (size_t i = 0; i < nc; i++)
				if (ndim == 2)
					Slice  (tmp, i, ft ->* Slice  (tmp, i));
				else
					Volume (tmp, i, ft ->* Volume (tmp, i));

#pragma omp for schedule (guided)
			
			// Antialiasing
			for (size_t x = 0; x < dims[0]; x++)
				for (size_t y = 0; y < dims[1]; y++)
					for (size_t z = 0; z < dims[2]; z++) {
						
						for (size_t c = 0; c < nc; c++) {
							
							ra [c] = (dims[2]-1) ? tmp (x, y, z, c) : tmp (x, y, c);
							
							size_t i = 0;
							for (size_t zi = 0; zi < af[2]; zi++)
								for (size_t yi = 0; yi < af[1]; yi++)
									for (size_t xi = 0; xi < af[0]; xi++, i++) {
									s (c, i) = (dims[2]-1) ?
											sens (x + xi * dims[0], y + yi * dims[1], z + zi * dims[2], c):
											sens (x + xi * dims[0], y + yi * dims[1],                   c);
									}
						}

						si = gemm (s,   s, 'C', 'N') ;

						if (treg > 0.0)
							si += reg;

						if (compgfm)
							gf = diag (si);

						si = inv (si);

						if (compgfm)
							gf = diag (si) * gf;

						si = gemm (si,  s, 'N', 'C');
						rp = gemm (si, ra, 'N', 'N');

						size_t i = 0;
						for (size_t zi = 0; zi < af[2]; zi++)
							for (size_t yi = 0; yi < af[1]; yi++)
								for (size_t xi = 0; xi < af[0]; xi++, i++) {
									res (x+xi*dims[0], y+yi*dims[1] , z+zi*dims[2]) = rp [i];
									//if (compgfm)
										//res (x, y + dims[1] * i, z, 1) = sqrt(abs(gf [i]));
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
	Matrix<CT>
	Trafo             (const Matrix<CT>& m) const {
		
		Matrix <CT> res;
		return res;


	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	operator* (const Matrix<CT>& m) const {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	operator->* (const Matrix<CT>& m) const {
		return Adjoint (m);
	}

private:

	/**
	 * @brief       Setup DFT operators
	 *
	 * @param  params Parameters
	 */
	inline void
	AllocateDFTs (const Params& params) {

		Workspace& w = Workspace::Instance();

		nthreads = params.Get<unsigned short>("nthreads");

		if (nthreads == 1)
			#pragma omp parallel default (shared)
			{
				nthreads = omp_get_num_threads ();
			}

		printf ("  allocating %zu %zu-dim ffts\n", nthreads, ndim);

		if (params.exists("nthreads"))
			nthreads = params.Get<unsigned short>("nthreads");

		m_dft = std::vector<DFT<T>*> (nthreads,(DFT<T>*)0);

		for (size_t i = 0; i < nthreads; i++)
			m_dft[i] = new DFT<T> (dims/*,
				w.Get<T>(params.Get<std::string>("mask_name")),
				w.Get<T>(params.Get<std::string>(  "pc_name")),
				w.Get<T>(params.Get<std::string>(  "b0_name"))*/);

	}


	/**
	 * @brief       Setup Tikhonov regularisation matrix
	 *
	 * @param  params Parameters
	 */
	inline void
	TikhonovMat (const Params& params) {

		treg = (params.exists("lambda")) ? params.Get<T>("lambda"): 0.0;
		printf ("  Tikhonov lambda (%.2e)\n", treg);
		assert (treg >= 0.0);
		if (treg > 0.0)
			reg = treg * eye<T>(aaf);

	}


	std::vector<DFT<T>*> m_dft;
	Matrix <CT>          sens;
	Matrix <size_t>      dims;
	int                  nthreads;
	Matrix <size_t>      af;
	size_t               ndim;
	size_t               nc;
	bool                 compgfm;
	T            		 treg;
	Matrix<T>            reg;
	bool                 initialised;
	size_t               aaf;

};


#endif
