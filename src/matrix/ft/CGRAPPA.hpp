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
#include "Workspace.hpp"
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
		m_af  = 1;
		m_nc  = 1;
	}


	/**
	 * @brief Construct with parameters
	 */
	CGRAPPA (const Params& p) : m_nc(1) {

		Workspace& ws = Workspace::Instance();
		printf ("Initialising GRAPPA operator ...\n");

		// Scan and ACS sizes
		CheckScan (p);
		CheckACS (p,ws);
		int err = 0;
		assert (m_sdims[0] == m_adims[0]);

		// Kernel dims usually ROxPExSS = 5x4x4
		CheckKernel (p);

		// Acceleartion factors
		m_af    = p.Get<Matrix<size_t> >("acc_factors");
		size_t af = prod(m_af);
		std::cout << err++ << std::endl;

		// Source and target matrices
		size_t nr = (m_adims[1]-m_kdims[1]-1) * (m_adims[2]-m_kdims[2]-1) * af;
		std::cout << err++ << std::endl;

		Matrix<CT> src = zeros<CT>(m_nc * prod(m_kdims), nr);
		std::cout << err++ << std::endl;
		Matrix<CT> trg = zeros<CT>(m_nc * af           , nr);

		std::cout << err++ << std::endl;

		std::cout << size(src);
		std::cout << err++ << std::endl;

		std::cout << size(trg);
		std::cout << err++ << std::endl;

		for (size_t y = floor(m_kdims[1]/2); y < m_adims[1]-floor(m_kdims[1]/2); ++y)
			for (size_t x = floor(m_kdims[0]/2); x < m_adims[0]-floor(m_kdims[0]/2); ++x) // f.e. 2-253
				for (size_t c = 0, /*srci = 0,*/ trgi = 0; c < m_nc; ++c) {
					for (size_t yt = y; yt < af; ++yt, ++trgi) {
						printf ("%zu, %zu, %zu, %zu, %zu\n", c, trgi, c, x, yt);
						trg (c,trgi) = m_acs (c,x,yt);
					}
					/*
					for (size_t ys = 0; ys < y+(m_kdims[1]-1)*af; ys += m_af[1]) {
						for (size_t xs = x-floor(m_kdims[0]/2); xs < x+floor(m_kdims[0]/2); ++xs, ++srci)
							src (c,srci) = m_acs (c,xs,ys);
					}*/
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


	bool CheckScan (const Params& p) {
		m_sdims = p.Get<Matrix<size_t> >("scan_dims");
		std::cout << "  Data dims: " << m_sdims;
		return true;
	}

	bool CheckKernel (const Params& p) {
		m_kdims = p.Get<Matrix<size_t> >("kern_dims");
		std::cout << "  Kernel dims: " << m_kdims;
		return true;
	}

	bool CheckACS (const Params& p, Workspace& ws) {
		bool integrated_acs = (p.exists("acs_integrated")) ? p.Get<bool>("integrated_acs") : false;
		assert((p.exists("acs_dims") && integrated_acs) || p.exists("acs_name"));
		m_acs   = integrated_acs ?
				  zeros<CT>(p.Get<Matrix<size_t> >("acs_dims")) :
				  ws.Get<CT>(p.Get<std::string>("acs_name"));
		m_adims = size(m_acs);
		m_nc    = m_adims[0];
		std::cout << "  ACS dims: " << m_adims;
		std::cout << "  NC: " << m_nc << std::endl;

		return true;
	}

	Matrix<CT>           m_weights; /**< @brief Correction patch  */
	Matrix<CT>           m_acs;     /**< @brief ACS lines         */

	Matrix<size_t>       m_kdims;   /**< @brief Kernel dimensions */
	Matrix<size_t>       m_adims;
	Matrix<size_t>       m_d;       /**< @brief Dimensions        */
	Matrix<size_t>       m_af;      /**< @brief Acceleration factors */
	Matrix<size_t>       m_sdims;   /**< @brief Scan dimensions */

	std::vector<DFT<T> > m_dft;     /**< @brief DFT operator      */

	size_t               m_nc;      /**< @brief Number of receive channels */


};

#endif /* __CGRAPPA_HPP__ */
