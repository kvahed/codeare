/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __REMOTE_CONNECTOR_H__
#define __REMOTE_CONNECTOR_H__

#include "Configurable.hpp"
#include "Connector.hpp"
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

	class LocalConnector;
	
	/**
	 * @brief               Remotely connected reconstruction client 
	 */
	class RemoteConnector : public Configurable {
		
		
	public:
		
		
		/**
		 * @brief           Construct and initialise remote interface
		 */
		RemoteConnector         (const char* name, const char* tracelevel);


		/**
		 * @brief           Construct from local connector (i.e. implicit cast)
		 */
		RemoteConnector         (const RRClient::LocalConnector&) ;


		/**
		 * @brief           Construct from local connector (i.e. implicit cast)
		 */
		RemoteConnector         (const RRClient::LocalConnector*&) ;


		/**
		 * @brief           Construct from local connector (i.e. implicit cast)
		 */
		RemoteConnector         (RRClient::LocalConnector&) ;


		/**
		 * @brief           Construct from local connector (i.e. implicit cast)
		 */
		RemoteConnector         (RRClient::LocalConnector*&) ;
		
		
		/**
		 * @brief           Clean up and destroy ORB
		 */
		~RemoteConnector        ();
		
		
 		/**
		 * @brief           Request data procession on remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Process
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		virtual error_code              
		Process             (const char* name);
		

 		/**
		 * @brief           Prepare backend
		 *
		 * @see             RRStrategy::ReconStrategy::Prepare
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		virtual error_code              
		Prepare             (const char* name);
		

 		/**
		 * @brief           Initialise remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Init
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		virtual error_code              
		Init                (const char* name);
		

 		/**
		 * @brief           Finalise remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Finalise
		 * @param  name     Recon method
		 * @return          Error error
		 */ 
		virtual error_code              
		Finalise            (const char* name);
		

		/**
		 * @brief           Transmit measurement data to remote service
		 *
		 * @see             Workspace::SetMatrix
		 * @param  name     Name
		 * @param  m        Complex data
		 */
		template <class T> void 
		SetMatrix           (const std::string& name, const Matrix<T>& m) const;
		
		
		/**
		 * @brief           Retrieve manipulated data from remote service
		 *
		 * @see             Workspace::GetMatrix
		 * @param  name     Name
		 * @param  m        Receive storage
		 */
		template <class T> void 
		GetMatrix           (const std::string& name, Matrix<T>& m) const;
		
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< @brief Remote Recon interface               */
		CORBA::ORB_var      m_orb;        /**< @brief Orb                                  */
		
		/**
		 * @brief           Get size from dimensions (Needed internally)
		 *
		 * @param  dims     Dimension array from the CORBA types 
		 * @return          Size
		 */
		long
		GetSize             (const longs dims) const;
		
		
	};
	
}


#endif // __REMOTE_CONNECTOR_H__
