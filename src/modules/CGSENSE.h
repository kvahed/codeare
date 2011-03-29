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

#ifndef __CGSENSE_H__
#define __CGSENSE_H__

#include "ReconStrategy.h"
#include "nfft3util.h"
#include "nfft3.h"

typedef std::complex<float> raw;

using namespace RRServer;


namespace RRStrategy {
	
	/**
	 * @brief Non uniform FFT
	 */
	class CGSENSE : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		CGSENSE  ();
		
		/**
		 * @brief Default destructor
		 */
		virtual 
			~CGSENSE () {};
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
			Process ();
		
		
		
	private:
		
		int                m_iter;
		int                m_verbose;
		
		Matrix < raw >     m_sens;
		Matrix < raw >     m_temp;
		
		nfft_plan          my_plan;            /* plan for nfft  */
		
		int*               N;
		int*               Nk;
		
		double*            kmax;
		
	};
	
};
#endif /* __CGSENSE_H__ */
