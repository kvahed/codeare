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
#include "Connection.hpp"
#include "Matrix.hpp"

#include "RRSModule.hh"

#include <complex>
#include <vector>





using namespace RRSModule;

/**
 * @brief Remote recon client
 */
namespace RRClient {

	template<class T> struct RemoteTraits;

	template<> struct RemoteTraits<float> {
		typedef rlfl_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_rlfl(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_rlfl(name.c_str(), ct);
		}
	};
	template<> struct RemoteTraits<double> {
		typedef rldb_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_rldb(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_rldb(name.c_str(), ct);
		}
	};
	template<> struct RemoteTraits<cxfl> {
		typedef cxfl_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_cxfl(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_cxfl(name.c_str(), ct);
		}
	};
	template<> struct RemoteTraits<cxdb> {
		typedef cxdb_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_cxdb(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_cxdb(name.c_str(), ct);
		}
	};
	template<> struct RemoteTraits<short> {
		typedef shrt_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_shrt(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_shrt(name.c_str(), ct);
		}
	};
	template<> struct RemoteTraits<long> {
		typedef long_data CORBA_Type;
		inline static void Send (const RRSInterface_var& m_rrsi, const std::string& name, const CORBA_Type& ct) {
			m_rrsi->set_long(name.c_str(), ct);
		}
		inline static void Retrieve (const RRSInterface_var& m_rrsi, const std::string& name, CORBA_Type& ct) {
			m_rrsi->get_long(name.c_str(), ct);
		}
	};


	/**
	 * @brief               Remotely connected reconstruction client 
	 */
	class RemoteConnector : public Configurable, public Connection {
		
		
	public:
		
		
		/**
		 * @brief           Construct and initialise remote interface
		 */
		RemoteConnector     (int argc, char** argv, const std::string& service_id,
				             const std::string& trace_level = "0", const std::string& client_id = "");


		/**
		 * @brief           Clean up and destroy ORB
		 */
		~RemoteConnector    ();
		
		
 		/**
		 * @brief           Request data procession on remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Process
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		codeare::error_code
		Process             (const std::string& name);
		

 		/**
		 * @brief           Prepare backend
		 *
		 * @see             RRStrategy::ReconStrategy::Prepare
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		codeare::error_code
		Prepare             (const std::string& name);
		

 		/**
		 * @brief           Initialise remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Init
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		codeare::error_code
		Init                (const std::string& name, const std::string& configuration);
		

 		/**
		 * @brief           Finalise remote service
		 *
		 * @see             RRStrategy::ReconStrategy::Finalise
		 * @param  name     Recon method
		 * @return          Error error
		 */ 
		codeare::error_code
		Finalise            (const std::string& name);
		

		/**
		 * @brief           Transmit measurement data to remote service
		 *
		 * @see             Workspace::SetMatrix
		 * @param  name     Name
		 * @param  m        Complex data
		 */
		template <class T> void 
		SetMatrix           (const std::string& name, const Matrix<T>& m) const {

			typename RemoteTraits<T>::CORBA_Type ct;

			size_t ms = m.Size();
			size_t nd = m.NDim();
			ct.dims.length(nd);
			ct.res.length(nd);

			for (int j = 0; j < nd; j++) {
				ct.dims[j] = m.Dim(j);
				ct.res[j]  = m.Res(j);
			}

			ct.vals.length(TypeTraits<T>::IsComplex() ? ms * 2 : ms);
			memcpy (&ct.vals[0], m.Ptr(), ms * sizeof(T));

			RemoteTraits<T>::Send(m_rrsi, name, ct);

		}

		
		
		/**
		 * @brief           Retrieve manipulated data from remote service
		 *
		 * @see             Workspace::GetMatrix
		 * @param  name     Name
		 * @param  m        Receive storage
		 */
		template <class T> void 
		GetMatrix           (const std::string& name, Matrix<T>& m) const {

			typename RemoteTraits<T>::CORBA_Type ct;

			RemoteTraits<T>::Retrieve (m_rrsi, name, ct);
			size_t nd = ct.dims.length();

			std::vector<size_t> mdims (nd);
			std::vector<float>  mress (nd);

			for (int j = 0; j < nd; j++) {
				mdims[j] = ct.dims[j];
				mress[j] = ct.res[j];
			}

			m = Matrix<T> (mdims, mress);
			size_t ms = m.Size();
			memcpy (&m[0], &ct.vals[0], ms * sizeof(T));

		}
		
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< @brief Remote Recon interface               */
		CORBA::ORB_var      m_orb;        /**< @brief Orb                                  */
        std::string         m_client_id;
		
		/**
		 * @brief           Get size from dimensions (Needed internally)
		 *
		 * @param  dims     Dimension array from the CORBA types 
		 * @return          Size
		 */
		long
		GetSize             (const longs& dims) const;
		
		
	};
	
}


#endif // __REMOTE_CONNECTOR_H__
