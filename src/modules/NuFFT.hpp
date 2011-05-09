/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#ifndef __NUFFT_HPP__
#define __NUFFT_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;

namespace RRStrategy {
	
	/**
	 * @brief Non uniform FFT
	 */
	class NuFFT : public ReconStrategy {
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		NuFFT  ();
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~NuFFT () {};
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Process ();
		
	private:
		
		int N [2];
		int Nk[2];
		
		// We only support 2 dims for the time being
		int m_dim;
		
		// Some more stuff
		int       m_N;
		int       m_iter;

		nfft_plan           m_fplan;            /**< nfft plan */
		solver_plan_complex m_iplan;

		
	};
	
}
#endif /* __NUFFT_HPP__ */
