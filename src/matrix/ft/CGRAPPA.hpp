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
			try {
				m_kernel = zeros<CT>(p.Get<Vector<size_t> >("kernel_size"));
			} catch (const std::exception& e) {
				std::cerr << "  WARNING - CGRAPPA: invalid kernel size definition, defaulting to 4x5" << std::endl;
			}
		} else {
			std::cerr << "  WARNING - CGRAPPA: kernel size unspecified, defaulting to 4x5" << std::endl;
			m_kernel = zeros<CT>(4,5);
		}

		// AC data
		if (p.exists("ac_data")) {
			try {
				m_ac_data = p.Get<Matrix<CT> >("ac_data");
				m_nc = size(m_ac_data,2); // # Coils
			} catch (const std::exception& e) {
				std::cerr << "  ERROR - CGRAPPA: auto calibration data is mandatory input.";
				assert (false);
			}
		}

		// Tikhonov lambda
		if (p.exists("lambda"))
			m_lambda = fp_cast(p["lambda"]);
		else
			m_lambda = T(0.);

	}
/*
	[sx,sy,nc] = size(data);

	tmp = im2row(data,kSize);
	[tsx,tsy,tsz] = size(tmp);
	A = reshape(tmp,tsx,tsy*tsz);

	AtA = A'*A;

	kernel = AtA;
	kernel = reshape(kernel,kSize(1),kSize(2),nc,size(kernel,2));
*/
	inline void CalcCalibMatrix () {
		Matrix<CT> tmp = Im2Row();
	}

	/*
	 *
	 function res = im2row(im, winSize)
%res = im2row(im, winSize)
[sx,sy,sz] = size(im);

res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
    end
end
	 *
	 */

	inline Matrix<CT> Im2Row () {
		Vector<size_t> ac_size = size(m_ac_data);
		Vector<size_t> kernel_size = size(m_kernel);
		Vector<size_t> calib_mat_size (3);
		calib_mat_size[0] = (ac_size[0] - kernel_size[0] + 1) * (ac_size[1] - kernel_size[1] + 1);
		calib_mat_size[1] = prod (kernel_size);
		calib_mat_size[2] = ac_size[2];
		Matrix<CT> ret (calib_mat_size);
		for (size_t j = 0, count = 0; j < kernel_size[1]; ++j)
			for (size_t i = 0; i < kernel_size[0]; ++i)
				for (size_t l = 0; l < calib_mat_size[0]; ++l)
					for (size_t m = 0; l < calib_mat_size[1]; ++m)
						for (size_t n = 0; n < ac_size[2]; ++n)
							ret (l, count, n) = m_ac_data(i+l,j+m,n);
		return ret;
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

	Matrix<CT>           m_weights; /**< @brief Correction patch     */
	Matrix<CT>           m_ac_data;     /**< @brief ACS lines            */
	Matrix<CT>           m_kernel;  /**< @brief GRAPPA kernel        */

	Matrix<size_t>       m_kdims;   /**< @brief Kernel dimensions    */
	Matrix<size_t>       m_adims;
	Matrix<size_t>       m_d;       /**< @brief Dimensions           */
	Matrix<size_t>       m_af;      /**< @brief Acceleration factors */
	Matrix<size_t>       m_sdims;   /**< @brief Scan dimensions      */

	T m_lambda;

	std::vector<DFT<T> > m_dft;     /**< @brief DFT operator         */

	size_t               m_nc;      /**< @brief Number of receive channels */

};

#endif /* __CGRAPPA_HPP__ */
