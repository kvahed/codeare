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
		 * @brief     Default constructor
		 */
		ReconServant  ();
		
		/**
		 * @brief     Default destructor
		 */
		virtual 
		~ReconServant ();
		
		/**
		 * @brief     Process startegy
		 *
		 * @param s   Initialised strategy
		 * @return    Sucess
		 */
		virtual error_code
		Process       (const short s);
		
		/**
		 * @brief     Initialise strategy
		 * 
		 * @param name Name of processing library
		 * @return    Id of 
		 */
		virtual short int
		Init          (const char* name);
		
		/**
		 * @brief     Finalise algorithm
		 */
		virtual error_code
		Finalise      (const short s);
		
		/**
		 * @brief     Get data from recon
		 */
		raw_data* 
		raw           ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		raw           (const raw_data&);
		
		/**
		 * @brief     Get data from recon
		 */
		raw_data* 
		rhelper           ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		rhelper           (const raw_data&);
		
		/**
		 * @brief     Get data from recon
		 */
		helper_data* 
		helper        ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		helper        (const helper_data&);
		
		/**
		 * @brief     Get data from recon
		 */
		helper_data* 
		kspace        ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		kspace        (const helper_data&);
		
		/**
		 * @brief     Get data from recon
		 */
		pixel_data* 
		pixel         ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		pixel         (const pixel_data&);
		
		/**
		 * @brief     Get data from recon
		 */
		char*
		config        ();
		
		/**
		 * @brief     Set data for recon
		 */
		void 
		config        (const char*);
		
		
	private:
		
		char*                 m_config;   /**< Configuration xml document  */
		std::vector<ReconContext*> m_contexts; /**< Reconstruction contexts     */

		
	};

}

#endif // __RECON_SERVANT_HPP__
