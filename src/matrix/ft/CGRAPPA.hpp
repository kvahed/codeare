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
	 * @brief Construct with parameters
	 */
	CGRAPPA (const Params& p) {

		// Number of coils
		if (p.exists("n_coils"))
			m_nc = unsigned_cast(p["n_coils"]);
		else {
			std::cerr << "  ERROR - CGRAPPA: Operator coil dimension definition" << std::endl;
			assert (false);
		}

		// Kernel size
		if (p.exists("kernel_size")) {
			/*try {
				m_kernel_size = p.Get<Matrix<size_t>&>("kernel_size");
			} catch (const std::exception& e) {
				std::cerr << "  ERROR - CGRAPPA: Operator needs kernel size definition" << std::endl;
				assert (false);
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
