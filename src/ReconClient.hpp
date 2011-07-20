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

#include "Matrix.hpp"
#include "Configurable.hpp"


#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

using namespace RRSModule;

/**
 * @brief Remote recon client
 */
namespace RRClient {
	
	/**
	 * @brief              CORBA ICE reconstruction client 
	 */
	class ReconClient :
	public Configurable {
		
		
	public:
		
		/**
		 * @brief           Construct and initialise remote recon interface
		 */
		ReconClient         (const char* name, const char* tracelevel);
		
		
		/**
		 * @brief           Default destructor
		 */
		~ReconClient        ();
		
		/**
		 * @brief           Default destructor
		 */
		void 
		Cleanup             ();
		
 		/**
		 * @brief           Request data procession on recon service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Process             (const char* name);
		
 		/**
		 * @brief           Request data procession on recon service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Init                (const char* name);
		
 		/**
		 * @brief           Request data procession on recon service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Finalise            (const char* name);
		
		/**
		 * @brief           Set raw my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		SetRaw              (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Return raw data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		GetRaw              (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Set raw my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		SetRHelper          (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Return raw data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		GetRHelper          (Matrix< std::complex<float> >& M);
		
		
		/**
		 * @brief           Set helper my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void
		SetHelper           (Matrix< double >& M);
		
		
		/**
		 * @brief           Return helper data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void
		GetHelper           (Matrix< double >& M);
		
		
		/**
		 * @brief           Set kspace my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void
		SetKSpace           (Matrix< double >& M);
		
		
		/**
		 * @brief           Return kspace data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void
		GetKSpace           (Matrix< double >& M);
		
		
		/**
		 * @brief           Set my Pixel data
		 * 
		 * @param  M        Given matrix
		 */
		void 
		SetPixel            (Matrix<short>& M);
		
		
		/**
		 * @brief           Put my Pixel data into 
		 * 
		 * @param  M        Given matrix
		 */
		void 
		GetPixel            (Matrix<short>& M);
		
		
		/**
		 * @brief           Get size from dimensions
		 *
		 * @return          Size
		 */
		long
		GetSize            (longs dims);
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< Remote Recon interface               */
		CORBA::ORB_var      m_orb;        /**< Orb                                  */
		std::vector<short>  m_rstrats;    /**< Remote reconstruction strategies    */
		
	};
	
}


#endif // __RECON_CLIENT_H__
