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

#ifndef __CGRAPPA_HPP__
#define __CGRAPPA_HPP__

#include "FT.hpp"

/**
 * @brief GRAPPA operator<br/>
 *        Griswold et al. MRM 2002, vol. 47 (6) pp. 1202-1210
 */
template <class T>
class CGRAPPA : public FT<T> {



public:


	/**
	 * @brief          Default constructor
	 */
	CGRAPPA           () {
		m_dft = 0;
		m_af  = 1;
		m_nc  = 1;
	}


	/**
	 * @brief          Construct with
	 *
	 * @param   sl     Side length of full images
	 * @param   nc     # receive channels
	 * @param   ksize  GRAPPA kernel size
	 * @param   nacs   # ACS lines
	 * @param   af     Acceleration factors
	 */
	CGRAPPA           (const Matrix< size_t >& sl, const Matrix<size_t>& kdims,
				       const size_t& nc, const size_t& nacs, const Matrix<size_t>& af) {

		m_dft    = 0;
		m_af     = af;
		m_nc     = nc;
		m_kdims  = kdims;

		m_dft    = new DFT<T>*[1];
		m_dft[0] = new DFT<T> (sl);

	}


	/**
	 * @brief    Clean up and destroy
	 */
	virtual
	~CGRAPPA () {};


	/**
	 * @brief    Adjoint transform
	 */
	Matrix< std::complex<T> >
	Adjoint (Matrix< std::complex<T> >) const {

		Matrix<T> res;
		return res;

	}


	/**
	 * @brief    Forward transform
	 */
	Matrix< std::complex<T> >
	Trafo (Matrix< std::complex<T> >) const {

		Matrix<T> res;
		return res;

	}

private:


	Matrix< std::complex<T> > m_weights; /**< @brief Correction patch  */
	Matrix< std::complex<T> > m_acs;     /**< @brief ACS lines         */

	Matrix< size_t >          m_kdims;   /**< @brief Kernel dimensions */
	Matrix< size_t >          m_d;       /**< @brief Dimensions        */
	Matrix< size_t >          m_af;      /**< @brief Acceleration factors */

	DFT<T>**                  m_dft;     /**< @brief DFT operator      */

	size_t                    m_nc;      /**< @brief Number of receive channels */


	/**
	 * @brief Compute GRAPPA weights
	 */
	void ComputeWeights () {

		ticks       tic     = getticks();
		printf ("  (Re-)Computing weights \n");

		// # Kernel repititios in ACS
		int nr = (acs_dim[1] - (kern_dim[1]-1) * R[1]) * (acs_dim[0] - (kern_dim[0] - 1));
		int ns = nc * kern_dim[0] * kern_dim[1]; // # Source points
		int nt = nc * (R[1]-1);      // # Target points

		Matrix<cxfl> s (ns, nr);  // Source
		Matrix<cxfl> t (nt, nr);  // Target

		int c = 0;


		//Grappa pattern for standard kernel 5x4
		Matrix<double> p = zeros<double> ((kern_dim[1]-1)*R[1]+1, kern_dim[0], nc);
		printf ("  patch size in ACS: %s\n", DimsToCString(p));

		// Source points
		for (int lin = 0; lin < kern_dim[1]; lin++)
			for (int col = 0; col < kern_dim[0]; col++)
				p (col,lin*R[1],0) = 1.0;

		// Target points
		size_t lins = 1 + (ceil(kern_dim[1]/2)-1) * R[1];
		for (int lin = lins; lin < lins + R[1] -1; lin++)
			p(ceil(p.Dim(0)/2), lin) = -1.0;


		printf ("  source matrix size %s\n", DimsToCString(s));

		int ni  = acs_dim[0] - d[0];
		int nj  = acs_dim[1] - (kern_dim[1]-1) * R[1];

		for (int i = d[0]; i < ni; i++)
			for (int j = 0, pos = 0; j < nj; j++, c++)
				for (int col = 0; col < kern_dim[0]; col++)
					for (int lin = 0; lin < kern_dim[1]; lin++)
						for (int ch = 0; ch < nc; ch++, pos++)
							s.At(pos + c*s.Dim(0)) = acs(i+col,j+lin,ch);

		//s.MXDump ("s.mat", "s", "");
		//p.MXDump ("p.mat", "p", "");
		//acs->MXDump ("acs.mat", "acs", "");


		//s = s.Pinv();

		/*
		  weights = t->*s;
		*/
		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	}


};

#endif /* __CGRAPPA_HPP__ */
