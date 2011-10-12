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

#ifndef __RECON_SERVANT_HPP__
#define __RECON_SERVANT_HPP__

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

#include <string>
#include <vector>

using namespace RRSModule;

/**
 * @brief Remote recon server
 */
namespace RRServer {
	
	class ReconContext;


	/**
	 * @brief Servant implementation 
	 *        Perform duties on remote server
	 */
	class ReconServant : 
		public POA_RRSModule::RRSInterface, 
		public PortableServer::RefCountServantBase {
		

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
		 * @param s    sth Initialised strategy
		 * @return     Sucess
		 */
		virtual error_code
		Process        (const short s);

		
		/**
		 * @brief      Prepare startegy (Needs initialisation @see Init)
		 *
		 * @param s    sth Initialised strategy
		 * @return     Sucess
		 */
		virtual error_code
		Prepare        (const short s);

		
		/**
		 * @brief      Initialise strategy (Configuration document needs to be set first @see config)
		 * 
		 * @param name Name of processing library
		 * @return     Id of 
		 */
		virtual short int
		Init          (const char* name);
		

		/**
		 * @brief     Finalise algorithm
		 *
		 * @param s   sth Intialised startegy
		 */
		virtual error_code
		Finalise      (const short s);
		

		/**
		 * @brief     Retreive measurement data
		 *
		 * @param  name  Name
		 * @return    Pointer to data
		 */
		void
		get_cplx      (const char* name, cplx_data&);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param  name  Name
		 * @param c   Complex measurement data
		 */
		void 
		set_cplx      (const char* name, const cplx_data& c);

		
		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param  name  Name
		 * @return    Pointer to real data
		 */
		void
		get_real      (const char* name, real_data& r);
		

		/**
		 * @brief     Set real data
		 *
		 * @param  name  Name
		 * @param r   Real helper data
		 */
		void 
		set_real      (const char* name, const real_data& r);
		

		/**
		 * @brief     Get pixel data from recon
		 *
		 * @param  name  Name
		 * @return    Pointer to pixel data
		 */
		void
		get_pixel     (const char* name, pixel_data& p);
		

		/**
		 * @brief     Set pixel data for recon
		 *
		 * @param  name  Name
		 * @param p   Pixel data
		 */
		void 
		set_pixel     (const char* name, const pixel_data& p);
		

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

		
		
	private:
		
		char*                      m_config;   /**< Serialised XML document  */
		std::vector<ReconContext*> m_contexts; /**< Reconstruction contexts (Abstraction layer to algorithms)*/
		
	};

}

#endif // __RECON_SERVANT_HPP__
