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

#include "Configurable.hpp"
#include "Matrix.hpp"


#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

#include <complex>
#include <vector>


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
		 * @brief           Prepare backend
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Prepare             (const char* name);
		

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
		 * @return          Error error
		 */ 
		error_code              
		Finalise            (const char* name);
		

		/**
		 * @brief           Transmit measurement data to remote service
		 *
		 * @param  name     Name
		 * @param  m        Complex data
		 */
		void 
		SetCplx             (const std::string name, Matrix<cplx>& m);
		
		
		/**
		 * @brief           Retrieve manipulated data from remote service
		 *
		 * @param  name     Name
		 * @param  m        Receive storage
		 */
		void 
		GetCplx             (const std::string name, Matrix<cplx>& m);
		
		
		/**
		 * @brief           Transmit real data to service (f.e. k-space weights)
		 *
		 * @param  name     Name
		 * @param  m        Real data
		 */
		void
		SetReal             (const std::string name, Matrix< double >& m);
		
		
		/**
		 * @brief           Get helper data after recon to ...
		 *
		 * @param  name     Name
		 * @param  m        Receive storage
		 */
		void
		GetReal             (const std::string name, Matrix< double >& m);
		
		
		/**
		 * @brief           Transmit gray-scale images to remote service
		 * 
		 * @param  name     Name
		 * @param  m        Short int data
		 */
		void 
		SetPixel            (const std::string name, Matrix<short>& m);
		
		
		/**
		 * @brief           Get gray-scale data after recon
		 * 
		 * @param  name     Name
		 * @param  m        Short int storage
		 */
		void 
		GetPixel            (const std::string name, Matrix<short>& m);
		
		
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
