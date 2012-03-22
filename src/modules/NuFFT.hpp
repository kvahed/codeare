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

#ifndef __NUFFT_HPP__
#define __NUFFT_HPP__

#include "NFFT.hpp"
#include "nfftstub.h"

#include "ReconStrategy.hpp"

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
		~NuFFT ();
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Process ();
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Prepare ();
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Init ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Finalise ();
		
	private:
		
		NFFT<cxdb> *nfft;

	};
	
}
#endif /* __NUFFT_HPP__ */
