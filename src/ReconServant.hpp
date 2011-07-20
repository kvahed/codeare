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

#include "Matrix.hpp"
#include "ReconContext.hpp"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif


using namespace RRSModule;

/**
 * @brief Remote recon server
 */
namespace RRServer {
	

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
		raw_data* 
		raw           ();
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param r   Complex measurement data
		 */
		void 
		raw           (const raw_data& r);

		
		/**
		 * @brief     Get complex helper data from recon
		 *
		 * @return    Pointer to data
		 */
		raw_data* 
		rhelper       ();
		

		/**
		 * @brief     Set complex helper data for recon
		 *
		 * @param r   Complex helper data (f.e. sensitivities)
		 */
		void 
		rhelper       (const raw_data& r);
		

		/**
		 * @brief     Get real helper data from recon
		 *
		 * @return    Pointer to real data
		 */
		helper_data* 
		helper        ();
		

		/**
		 * @brief     Set real helper data for recon
		 *
		 * @param r   Real helper data
		 */
		void 
		helper        (const helper_data& r);
		

		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param     Pointer to real helper data
		 */
		helper_data* 
		kspace        ();

		
		/**
		 * @brief     Set kspace for recon
		 *
		 * @param k   Kspace 
		 */
		void 
		kspace        (const helper_data& k);

		
		/**
		 * @brief     Get pixel data from recon
		 *
		 * @return    Pointer to pixel data
		 */
		pixel_data* 
		pixel         ();
		

		/**
		 * @brief     Set pixel data for recon
		 *
		 * @param p   Pixel data
		 */
		void 
		pixel         (const pixel_data& p);
		

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
		 * @param     Configuration
		 */
		void 
		config        (const char*);

		
		
	private:
		
		char*                      m_config;   /**< Serialised XML document  */
		std::vector<ReconContext*> m_contexts; /**< Reconstruction contexts (Abstraction layer to algorithms)*/
		
	};

}

#endif // __RECON_SERVANT_HPP__
