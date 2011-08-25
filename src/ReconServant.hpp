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
		public POA_RRSModule::RRSInterface , 
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
		 * @return    Pointer to data
		 */
		cplx_data* 
		get_cplx      (const char* name);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param r   Complex measurement data
		 */
		void 
		set_cplx      (const char* name, const cplx_data& r);

		
		/**
		 * @brief     Get real helper data from recon
		 *
		 * @return    Pointer to real data
		 */
		real_data* 
		get_real      (const char* name);
		

		/**
		 * @brief     Set real data
		 *
		 * @param r   Real helper data
		 */
		void 
		set_real      (const char* name, const real_data& r);
		

		/**
		 * @brief     Get pixel data from recon
		 *
		 * @return    Pointer to pixel data
		 */
		pixel_data* 
		get_pixel     (const char* name);
		

		/**
		 * @brief     Set pixel data for recon
		 *
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
