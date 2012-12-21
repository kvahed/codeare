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

#ifndef __RECON_SERVANT_HPP__
#define __RECON_SERVANT_HPP__

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

#include <string>
#include <map>
#include "ReconContext.hpp"
#include "FunctorContainer.hpp"

using namespace RRSModule;
using namespace RRStrategy;
using namespace std;

/**
 * @brief Remote recon server
 */


namespace RRServer {
	
	/**
	 * @brief Servant implementation <br/>
	 *        Perform duties on remote server
	 */
	class ReconServant : 
		public POA_RRSModule::RRSInterface, 
		public PortableServer::RefCountServantBase,
		public FunctorContainer {
		
		
	public:
		
		
		/**
		 * @brief     Construct and prepare configuration document
		 */
		ReconServant  ();
		
		
		/**
		 * @brief     Default destructor
		 */
		virtual 
		~ReconServant ();
		
		
		/**
		 * @brief      Process startegy (Needs initialisation @see Init)
		 *
		 * @see        RRStrategy::ReconStrategy::Process()
		 *
		 * @param name Name of library
		 * @return     Sucess
		 */
		virtual error_code
		Process        (const char* name);
		
		
		/**
		 * @brief      Prepare startegy (Needs initialisation @see Init)
		 *
		 * @see        RRStrategy::ReconStrategy::Prepare()
		 *
		 * @param name Name of library
		 * @return     Sucess
		 */
		virtual error_code
		Prepare        (const char* name);
		
		
		/**
		 * @brief      Initialise strategy (Configuration document needs to be set first @see config)
		 * 
		 * @see        RRStrategy::ReconStrategy::Init()
		 *
		 * @param name Name of processing library
		 * @return     success
		 */
		virtual error_code
		Init          (const char* name);
		

		/**
		 * @brief     Finalise algorithm
		 *
		 * @see        RRStrategy::ReconStrategy::Finalise()
		 *
		 * @param name Name of library
		 */
		virtual error_code
		Finalise      (const char* name = 0);
		

		/**
		 * @brief     Clean up left over objects
		 *
		 * @return    Success
		 */
		virtual error_code
		CleanUp       ();
		

		/**
		 * @brief       Retreive measurement data
		 *
		 * @param  name Name
		 * @param  c    Data respository
		 */
		void
		get_cxfl      (const char* name, cxfl_data& c);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param  name  Name
		 * @param c   Complex measurement data
		 */
		void 
		set_cxfl      (const char* name, const cxfl_data& c);

		
		/**
		 * @brief       Retreive measurement data
		 *
		 * @param  name Name
		 * @param  c    Data respository
		 */
		void
		get_cxdb      (const char* name, cxdb_data& c);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param  name  Name
		 * @param c   Complex measurement data
		 */
		void 
		set_cxdb      (const char* name, const cxdb_data& c);

		
		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param  name  Name
		 * @param  r    Data
		 */
		void
		get_rldb      (const char* name, rldb_data& r);
		

		/**
		 * @brief     Set real data
		 *
		 * @param  name  Name
		 * @param r   Real helper data
		 */
		void 
		set_rldb      (const char* name, const rldb_data& r);
		

		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param  name  Name
		 * @param  r    Data
		 */
		void
		get_rlfl      (const char* name, rlfl_data& r);
		

		/**
		 * @brief     Set real data
		 *
		 * @param  name  Name
		 * @param r   Rlfl data
		 */
		void 
		set_rlfl      (const char* name, const rlfl_data& r);
		

		/**
		 * @brief     Get shrt data from recon
		 *
		 * @param  name  Name
		 * @param  p     Data
		 */
		void
		get_shrt     (const char* name, shrt_data& p);
		

		/**
		 * @brief     Set shrt data for recon
		 *
		 * @param  name  Name
		 * @param p   Shrt data
		 */
		void 
		set_shrt     (const char* name, const shrt_data& p);
		

		/**
		 * @brief     Get long data from recon
		 *
		 * @param  name  Name
		 * @param  p     Data
		 */
		void
		get_long     (const char* name, long_data& p);
		

		/**
		 * @brief     Set long data for recon
		 *
		 * @param  name  Name
		 * @param p   Long data
		 */
		void 
		set_long     (const char* name, const long_data& p);
		

		/**
		 * @brief     Get serialised configuration from backend
		 *
		 * @return    Serialised configuration
		 */
		char*
		config        ();

		
		/**
		 * @brief     Set serialised configuration
		 *
		 * @param c   Serialised Configuration
		 */
		void 
		config        (const char* c);

		
	};

}

#endif // __RECON_SERVANT_HPP__
