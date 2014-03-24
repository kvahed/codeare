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

#ifndef __LOCAL_CONNECTOR_H__
#define __LOCAL_CONNECTOR_H__

#include <complex>
#include <vector>

#include "Configurable.hpp"
#include "ReconContext.hpp"
#include "Connection.hpp"
#include "Queue.hpp"

using namespace RRStrategy;

/**
 * @brief Recon client
 */
namespace RRClient {
	
	/**
	 * @brief               Locally connected reconstruction client 
	 */
	class LocalConnector :
        public Configurable, public Connection, public Queue {
		
		typedef std::map<std::string, ReconContext*> context_map;

		
	public:
		
		/**
		 * @brief       Default constructor
		 */
		LocalConnector ();


		/**
		 * @brief       Default constructor
		 */
		LocalConnector (int args, char** argv, const char* name, const char* debug);


		/**
		 * @brief       Disconnect
		 */
		virtual ~LocalConnector     ();
		
		
		/**
		 * @brief      Process startegy (Needs initialisation @see Init)
		 *
		 * @param name Name of processing library
		 * @return     Sucess
		 */
		virtual short
		Process        (const char* name);

		
		/**
		 * @brief      Prepare startegy (Needs initialisation @see Init)
		 *
		 * @param name Name of processing library
		 * @return     Sucess
		 */
		virtual short
		Prepare        (const char* name);

		
		/**
		 * @brief      Initialise strategy (Configuration document needs to be set first @see config)
		 * 
		 * @param name Name of processing library
		 * @return     Success
		 */
		virtual short
		Init           (const char* name, const char* configuration);
		

		/**
		 * @brief      Finalise algorithm
		 *
		 * @param name Name of processing library
		 * @return     Success
		 */
		virtual short
		Finalise       (const char* name = 0);
		

		/**
		 * @brief      Clean up left over objects
		 *
		 * @return     Success
		 */
		virtual short
		CleanUp        ();
		

		/**
		 * @brief Transmit measurement data to remote service
		 *
		 * @param  name     Name
		 * @param  m        Data
		 */
		template <class T> void 
		SetMatrix           (const std::string& name, Matrix<T>& m) const {			
			Workspace::Instance().SetMatrix(name, m);
		}
		
		
		/**
		 * @brief           Retrieve manipulated data from remote service
		 *
		 * @param  name     Name
		 * @param  m        Receive storage
		 */
		template <class T> codeare::error_code
		GetMatrix           (const std::string& name, Matrix<T>& m) const {
			return Workspace::Instance().GetMatrix(name, m);
		}
		
		
		
	private:
		
		context_map m_contexts; /**< Reconstruction contexts (Abstraction layer to algorithms)*/
		
		
	};
	
}


#endif // __LOCAL_CONNECTOR_H__
