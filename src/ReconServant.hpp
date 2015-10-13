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

#include "RRSModule.hh"

#include <string>
#include <map>
#include "ReconContext.hpp"
#include "Queue.hpp"
#include "omniORB4/omniInterceptors.h"

using namespace std;

/**
 * @brief Remote recon server
 */


namespace RRServer {


	template<class T> struct RemoteTraits;
	
	template<> struct RemoteTraits<RRSModule::rlfl_data> {
		typedef float Type;
	};
	template<> struct RemoteTraits<RRSModule::rldb_data> {
		typedef double Type;
	};
	template<> struct RemoteTraits<RRSModule::cxfl_data> {
		typedef cxfl Type;
	};
	template<> struct RemoteTraits<RRSModule::cxdb_data> {
		typedef cxdb Type;
	};
	template<> struct RemoteTraits<RRSModule::shrt_data> {
		typedef short Type;
	};
	template<> struct RemoteTraits<RRSModule::long_data> {
		typedef long Type;
	};


	/**
	 * @brief Servant implementation <br/>
	 *        Perform duties on remote server
	 */
	class ReconServant : 
		public POA_RRSModule::RRSInterface, 
		public PortableServer::RefCountServantBase,
		public Queue {
		
		

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
		virtual short
		Process        (const char* name);
		
		/**
		 * @brief      Prepare startegy (Needs initialisation @see Init)
		 *
		 * @see        RRStrategy::ReconStrategy::Prepare()
		 *
		 * @param name Name of library
		 * @return     Sucess
		 */
		virtual short
		Prepare        (const char* name);
		
		
		/**
		 * @brief      Initialise strategy (Configuration document needs to be set first @see config)
		 * 
		 * @see        RRStrategy::ReconStrategy::Init()
		 *
		 * @param name Name of processing library
		 * @return     success
		 */
		virtual short
		Init          (const char* name, const char* config, const char* client_id);
		

		/**
		 * @brief     Finalise algorithm
		 *
		 * @see        RRStrategy::ReconStrategy::Finalise()
		 *
		 * @param name Name of library
		 */
		virtual short
		Finalise      (const char* name = 0);
		

		/**
		 * @brief     Clean up left over objects
		 *
		 * @return    Success
		 */
		virtual short
		CleanUp       ();


		template <class CORBA_Type> void
		SetMatrix (const char* name, const CORBA_Type& c) {

			typedef typename RemoteTraits<CORBA_Type>::Type T;

			size_t nd = c.dims.length();
			Vector<size_t> mdims (nd);
			Vector<float>  mress (nd);

			for (size_t i = 0; i < nd; i++) {
				mdims[i] = c.dims[i];
				mress[i] = c.res[i];
			}

			Matrix<T> pm (mdims, mress);
			memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(T));

			Workspace::Instance().SetMatrix(name, pm);

		}

		template <class CORBA_Type> void GetMatrix (const char* name, CORBA_Type& c) {

			typedef typename RemoteTraits<CORBA_Type>::Type T;

			Matrix<T> tmp = Workspace::Instance().Get<T> (name);
			size_t cpsz = tmp.Size();
			size_t nd = tmp.NDim();
			c.dims.length(nd);
			c.res.length (nd);

			for (size_t j = 0; j < nd; j++) {
				c.dims[j] = tmp.Dim(j);
				c.res[j]  = tmp.Res(j);
			}

			c.vals.length(TypeTraits<T>::IsComplex() ? 2 * cpsz : cpsz);
			memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(T));

		}



		/**
		 * @brief       Retreive measurement data
		 *
		 * @param  name Name
		 * @param  c    Data respository
		 */
		void
		get_cxfl      (const char* name, RRSModule::cxfl_data& c);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param  name  Name
		 * @param c   Complex measurement data
		 */
		void 
		set_cxfl      (const char* name, const RRSModule::cxfl_data& c);

		
		/**
		 * @brief       Retreive measurement data
		 *
		 * @param  name Name
		 * @param  c    Data respository
		 */
		void
		get_cxdb      (const char* name, RRSModule::cxdb_data& c);
		

		/**
		 * @brief     Set measurement data
		 *
		 * @param  name  Name
		 * @param c   Complex measurement data
		 */
		void 
		set_cxdb      (const char* name, const RRSModule::cxdb_data& c);

		
		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param  name  Name
		 * @param  r    Data
		 */
		void
		get_rldb      (const char* name, RRSModule::rldb_data& r);
		

		/**
		 * @brief     Set real data
		 *
		 * @param  name  Name
		 * @param r   Real helper data
		 */
		void 
		set_rldb      (const char* name, const RRSModule::rldb_data& r);
		

		/**
		 * @brief     Get real helper data from recon
		 *
		 * @param  name  Name
		 * @param  r    Data
		 */
		void
		get_rlfl      (const char* name, RRSModule::rlfl_data& r);
		

		/**
		 * @brief     Set real data
		 *
		 * @param  name  Name
		 * @param r   Rlfl data
		 */
		void 
		set_rlfl      (const char* name, const RRSModule::rlfl_data& r);
		

		/**
		 * @brief     Get shrt data from recon
		 *
		 * @param  name  Name
		 * @param  p     Data
		 */
		void
		get_shrt     (const char* name, RRSModule::shrt_data& p);
		

		/**
		 * @brief     Set shrt data for recon
		 *
		 * @param  name  Name
		 * @param p   Shrt data
		 */
		void 
		set_shrt     (const char* name, const RRSModule::shrt_data& p);
		

		/**
		 * @brief     Get long data from recon
		 *
		 * @param  name  Name
		 * @param  p     Data
		 */
		void
		get_long     (const char* name, RRSModule::long_data& p);
		

		/**
		 * @brief     Set long data for recon
		 *
		 * @param  name  Name
		 * @param p   Long data
		 */
		void 
		set_long     (const char* name, const RRSModule::long_data& p);
		

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



        static void
        inform        (omni::omniInterceptors::assignUpcallThread_T::info_T &info);

		
	};

}

#endif // __RECON_SERVANT_HPP__
