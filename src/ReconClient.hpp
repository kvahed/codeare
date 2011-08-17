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

#ifndef __RECON_CLIENT_H__
#define __RECON_CLIENT_H__

#include <complex>
#include <vector>
#include "Configurable.hpp"


#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif


template <class T> class Matrix;

using namespace RRSModule;

/**
 * @brief Remote recon client
 */
namespace RRClient {
	
	/**
	 * @brief               Remote reconstruction client 
	 */
	class ReconClient :
	public Configurable {
		
		
	public:
		
		/**
		 * @brief           Construct and initialise remote interface
		 */
		ReconClient         (const char* name, const char* tracelevel);
		
		
		/**
		 * @brief           Destroy ORB
		 */
		~ReconClient        ();
		
		
 		/**
		 * @brief           Request data procession on remote service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Process             (const char* name);
		

 		/**
		 * @brief           Initialise remote service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Init                (const char* name);
		

 		/**
		 * @brief           Finalise remote service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Finalise            (const char* name);
		

		/**
		 * @brief           Transmit measurement data to remote service
		 *
		 * @param  M        Complex data
		 */
		void 
		SetRaw              (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Retrieve manipulated data from remote service
		 *
		 * @param  M        Receive storage
		 */
		void 
		GetRaw              (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Transmit complex helper data (f.e. sensitivity maps)
		 *
		 * @param  M        Complex data
		 */
		void 
		SetRHelper          (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Get complex helper data (f.e. 2nd set of images)
		 *
		 * @param  M        Receive storage
		 */
		void 
		GetRHelper          (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Transmit real data to service (f.e. k-space weights)
		 *
		 * @param  M        Real data
		 */
		void
		SetHelper           (Matrix< double >& M);
		
		
		/**
		 * @brief           Get helper data after recon to ...
		 *
		 * @param  M        Receive storage
		 */
		void
		GetHelper           (Matrix< double >& M);
		
		
		/**
		 * @brief           Transmit k-space data to storage
		 *
		 * @param  M        Real data
		 */
		void
		SetKSpace           (Matrix< double >& M);
		
		
		/**
		 * @brief           Get k-space data after manipulation?
		 *
		 * @param  M        Real data
		 */
		void
		GetKSpace           (Matrix< double >& M);
		
		
		/**
		 * @brief           Transmit gray-scale images to remote service
		 * 
		 * @param  M        Short int data
		 */
		void 
		SetPixel            (Matrix<short>& M);
		
		
		/**
		 * @brief           Get gray-scale data after recon
		 * 
		 * @param  M        Short int storage
		 */
		void 
		GetPixel            (Matrix<short>& M);
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< Remote Recon interface               */
		CORBA::ORB_var      m_orb;        /**< Orb                                  */
		std::vector<short>  m_rstrats;    /**< Remote reconstruction strategies    */
		
		/**
		 * @brief           Get size from dimensions (Needed internally)
		 *
		 * @param  dims     Dimension array from the CORBA types 
		 * @return          Size
		 */
		long
		GetSize             (longs dims);
		
		
	};
	
}


#endif // __RECON_CLIENT_H__
