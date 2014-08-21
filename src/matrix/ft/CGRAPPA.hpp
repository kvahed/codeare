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
#include "Lapack.hpp"
#include "Print.hpp"
#include "Access.hpp"


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
	CGRAPPA() :  m_nthreads(1), m_lambda(0), m_nc(1) {}


	/**
	 * @brief Construct with parameters
	 */
	CGRAPPA (const Params& p) {

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
		std::cout << "  kernel size: " << size(m_kernel) << std::endl;

		// AC data
		if (p.exists("ac_data")) {
			try {
				m_ac_data = p.Get<Matrix<CT> >("ac_data");
			} catch (const std::exception& e) {
				std::cerr << "  ERROR - CGRAPPA: auto calibration data is mandatory input.";
				assert (false);
			}
		}
		m_nc = size(m_ac_data,2); // # Coils
		std::cout << "\n  # coils: " << m_nc << std::endl;

		// Tikhonov lambda
		if (p.exists("lambda"))
			m_lambda = fp_cast(p["lambda"]);
		else
			m_lambda = T(0.);

		// Parallelisation
		if (p.exists("nthreads"))
			m_nthreads = fp_cast(p["threads"]);
		else
			m_nthreads = 1;

		CalcCalibMatrix();

	}

//	printf( "ret(%zu, %zu, %zu, %zu) = m_ac_data(%zu, %zu, %zu)\n", l, m, count, n, i+l, j+m, n); fflush(stdout);



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
		Matrix<CT> res = kspace;
		//for (size_t coil = 0; coil < m_nc; ++coil)
			//Slice(res, coil, ReconstructSingle(Slice(kspace, coil)));

		return zpad(kspace, 300, 300, 10);
	}


	/**
	 * @brief    Forward transform
	 */
	Matrix<CT>
	Trafo (const Matrix<CT>& image) const {
		Matrix<CT> res;
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

	/*
	function [res] = ARC(kData, AtA, kSize, c,lambda)

	[sx,sy,nCoil] = size(kData);


	kData = zpad(kData,[sx+kSize(1)-1, sy+kSize(2)-1,nCoil]);
	dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK((end+1)/2,(end+1)/2,c) = 1;
	idxy = find(dummyK);
	res = zeros(sx,sy);

	MaxListLen = 100;
	LIST = zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
	KEY =  zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
	count = 0;

	%H = waitbar(0);
	for y = 1:sy
		for x=1:sx
		%	waitbar((x + (y-1)*sx)/sx/sy,H);
			tmp = kData(x:x+kSize(1)-1,y:y+kSize(2)-1,:);
			pat = abs(tmp)>0;
			if pat(idxy) | sum(pat)==0
				res(x,y) = tmp(idxy);
			else
				key = pat(:);
	            idx = 0;
				for nn=1:size(KEY,2);
					if sum(key==KEY(:,nn))==length(key)
					   idx = nn;
					   break;
				   	end
				end
				if idx == 0
					count = count + 1;
					kernel = calibrate(AtA,kSize,nCoil,c,lambda,pat);
					KEY(:,mod(count,MaxListLen)+1) = key(:);
					LIST(:,mod(count,MaxListLen)+1) = kernel(:);
					%disp('add another key');size(KEY,2)
				else
					kernel = LIST(:,idx);
				end
				res(x,y) = sum(kernel(:).*tmp(:));
			end

		end
	end
*/
	inline const Matrix<CT>& ReconstructSingle (const Matrix<CT>& coil) const {
		return coil;
	}

	inline void CalcCalibMatrix () {
		Vector<size_t> ac_size = size(m_ac_data);
		Vector<size_t> kernel_size = size(m_kernel);
		Vector<size_t> calib_mat_size (4);
		calib_mat_size[0] = ac_size[0] - kernel_size[0] + 1;
		calib_mat_size[1] = ac_size[1] - kernel_size[1] + 1;
		calib_mat_size[2] = prod (kernel_size);
		calib_mat_size[3] = ac_size[2];
		m_coil_calib = Matrix<CT> (calib_mat_size);
 		for (size_t j = 0, count = 0; j < kernel_size[1]; ++j)
			for (size_t i = 0; i < kernel_size[0]; ++i, ++count)
				for (size_t m = 0; m < calib_mat_size[1]; ++m)
					for (size_t l = 0; l < calib_mat_size[0]; ++l)
						for (size_t n = 0; n < ac_size[2]; ++n)
							m_coil_calib (l, m, count, n) = m_ac_data(i+l,j+m,n);
 		calib_mat_size[0] *= calib_mat_size[1];
 		calib_mat_size[1]  = calib_mat_size[2] * calib_mat_size[3];
 		calib_mat_size.PopBack();
 		calib_mat_size.PopBack();
 		m_coil_calib = resize(m_coil_calib, calib_mat_size);
		m_coil_calib = gemm (m_coil_calib, m_coil_calib, 'C');
	}

	Matrix<CT>           m_weights; /**< @brief Correction patch     */
	Matrix<CT>           m_ac_data;     /**< @brief ACS lines            */
	Matrix<CT>           m_kernel;  /**< @brief GRAPPA kernel        */
	Matrix<CT>           m_coil_calib;

	Matrix<size_t>       m_kdims;   /**< @brief Kernel dimensions    */
	Matrix<size_t>       m_adims;
	Matrix<size_t>       m_d;       /**< @brief Dimensions           */
	Matrix<size_t>       m_af;      /**< @brief Acceleration factors */
	Matrix<size_t>       m_sdims;   /**< @brief Scan dimensions      */

	T m_lambda;

	std::vector<DFT<T> > m_dft;     /**< @brief DFT operator         */

	size_t               m_nc;      /**< @brief Number of receive channels */
	size_t               m_nthreads;

};

#endif /* __CGRAPPA_HPP__ */
