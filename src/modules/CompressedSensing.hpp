/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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

#ifndef __COMPRESSED_SENSING_HPP__
#define __COMPRESSED_SENSING_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief CS reconstruction based on Sparse MRI v0.2 by Michael Lustig
	 */
	class CompressedSensing : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		CompressedSensing  () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~CompressedSensing () {};
		
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Process ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code 
		Init ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Finalise () {};


	private:
		
		int            m_dim;  /**< Image recon dim */
		int            m_N[3]; /**< Data side lengths */
		int            m_fft;  /**< FFT mode: 0 - Cartesian, 1 - Non Cartesian */
		int            m_maxiter; /**< Maximum # iterations */

		double         m_tvw;  /**< Total variation weights */
		double         m_xfmw; /**< Transform l1 penalty */

		//bwave::Wavelet m_wl;   /**< Wavelet */
		
		
	};

}
#endif /* __COMPRESSED_SENSING_H__ */
