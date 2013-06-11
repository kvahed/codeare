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

#include "DFT.hpp"
#include "GRAPPATraits.hpp"
#include "Print.hpp"


/**
 * @brief GRAPPA operator<br/>
 *        Griswold et al. MRM 2002, vol. 47 (6) pp. 1202-1210
 */
template <class T>
class CGRAPPA : public FT<T> {


	typedef std::complex<T> CT;


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
	 * @brief Construct with parameters
	 */
	CGRAPPA (const Params& p) : m_nc(1), m_dft(0) {

		Workspace& ws = Workspace::Instance();
		printf ("Initialising GRAPPA operator ...\n");

		// Scan and ACS sizes
		CheckScan(p, ws);
		CheckACS(p, ws);
		assert (m_sdims[0] == m_adims[0]);
		
		// Kernel dims usually ROxPExSS = 5x4x4
		m_kdims = p.Get<Matrix<size_t> >("kern_dims");

		// Acceleartion factors
		m_af    = p.Get<Matrix<size_t> >("acc_factors");
		size_t af = prod(m_af);

		// Source and target matrices
		size_t nr = (m_adims[1]-m_kdims[1]-1) * (m_adims[2]-m_kdims[2]-1) * af;
		Matrix<CT> src = zeros<CT>(m_nc * prod(m_kdims), nr);
		Matrix<CT> trg = zeros<CT>(m_nc * af           , nr);

		for (size_t y = floor(m_kdims[1]/2); y < m_adims[1]-floor(m_kdims[1]/2); ++y)
			for (size_t x = floor(m_kdims[0]/2); x < m_adims[0]-floor(m_kdims[0]/2); ++x) // f.e. 2-253
				for (size_t c = 0, srci = 0, trgi = 0; c < m_nc; ++c) {
					for (size_t yt = y; yt < af; ++yt, ++trgi)

						trg (c,trgi) = m_acs (c,x,yt);
					for (size_t ys = 0; ys < y+(m_kdims[1]-1)*af; ys += m_af[1])
						for (size_t xs = x-floor(m_kdims[0]/2); xs < x+floor(m_kdims[0]/2); ++xs, ++srci)
							src (c,srci) = m_acs (c,xs,ys);
				}



	}


	/**
	 * @brief    Clean up and destroy
	 */
	virtual
	~CGRAPPA () {};


	/**
	 * @brief    Adjoint transform
	 */
	Matrix<CT>
	Adjoint (const Matrix<CT>& kspace) const {

		Matrix<T> res;
		return res;

	}


	/**
	 * @brief    Forward transform
	 */
	Matrix<CT>
	Trafo (const Matrix<CT>& image) const {

		Matrix<T> res;
		return res;

	}

private:


	bool CheckScan (const Params& p, const Workspace& ws) {
		m_sdims = p.Get<Matrix<size_t> >("scan_dims");
		std::cout << "  Data dims: " << m_sdims;
		return true;
	}

	bool CheckACS (const Params& p, Workspace& ws) {
		bool integrated_acs = (p.exists("acs_integrated")) ? p.Get<bool>("integrated_acs") : false;
		assert((p.exists("acs_dims") && p.Get<bool>("acs_integrated")) || p.exists("acs_name"));
		m_acs   = (integrated_acs) ?
				  zeros<CT>(p.Get<Matrix<size_t> >("acs_dims")) :
				  ws.Get<CT>(p.Get<std::string>("acs_name"));
		m_adims = size(m_acs);
		m_nc    = m_adims[0];
		std::cout << "  ACS dims: " << size(m_acs);
		return true;
	}

	Matrix<CT>           m_weights; /**< @brief Correction patch  */
	Matrix<CT>           m_acs;     /**< @brief ACS lines         */

	Matrix<size_t>       m_kdims;   /**< @brief Kernel dimensions */
	Matrix<size_t>       m_adims;
	Matrix<size_t>       m_d;       /**< @brief Dimensions        */
	Matrix<size_t>       m_af;      /**< @brief Acceleration factors */
	Matrix<size_t>       m_sdims;   /**< @brief Scan dimensions */

	std::vector<DFT<T>*> m_dft;     /**< @brief DFT operator      */

	size_t               m_nc;      /**< @brief Number of receive channels */


	/**
	 * @brief Compute GRAPPA weights
	 */
	/*
	void ComputeWeights () {

		ticks       tic     = getticks();
		printf ("  (Re-)Computing weights \n");

		// # Kernel repititios in ACS
		int nr = (acs_dim[1] - (kern_dim[1]-1) * m_af[1]) * (acs_dim[0] - (kern_dim[0] - 1));
		int ns = m_nc * kern_dim[0] * kern_dim[1]; // # Source points
		int nt = m_nc * (m_af[1]-1);      // # Target points

		Matrix<cxfl> s (ns, nr);  // Source
		Matrix<cxfl> t (nt, nr);  // Target

		int c = 0;

		//Grappa pattern for standard kernel 5x4
		Matrix<double> p = zeros<double> ((kern_dim[1]-1)*m_af[1]+1, kern_dim[0], m_nc);
		printf ("  patch size in ACS: %s\n", DimsToCString(p));

		// Source points
		for (int lin = 0; lin < kern_dim[1]; lin++)
			for (int col = 0; col < kern_dim[0]; col++)
				p (col,lin*m_af[1],0) = 1.0;

		// Target points
		size_t lins = 1 + (ceil(kern_dim[1]/2)-1) * m_af[1];
		for (int lin = lins; lin < lins + m_af[1] -1; lin++)
			p(ceil(p.Dim(0)/2), lin) = -1.0;

		printf ("  source matrix size %s\n", DimsToCString(s));

		int ni  = acs_dim[0] - d[0];
		int nj  = acs_dim[1] - (kern_dim[1]-1) * m_af[1];

		for (int i = d[0]; i < ni; i++)
			for (int j = 0, pos = 0; j < nj; j++, c++)
				for (int col = 0; col < kern_dim[0]; col++)
					for (int lin = 0; lin < kern_dim[1]; lin++)
						for (int ch = 0; ch < nc; ch++, pos++)
							s.At(pos + c*s.Dim(0)) = acs(i+col,j+lin,ch);


		s = s.Pinv();
		weights = t->*s;


		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	}*/


};

#endif /* __CGRAPPA_HPP__ */
