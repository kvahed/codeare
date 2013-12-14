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

#include "common.h"
#include "RemoteConnector.hpp"
#include "CorbaExceptions.hpp"


#include <assert.h>
#include <complex>
#include <string>
#include <sstream>

namespace RRClient {


	RemoteConnector::RemoteConnector          (const std::string& service_id, const std::string& debug_level, const std::string& client_id) {
		
		try {

            if (!client_id.empty()) {
                m_client_id = client_id;
            } else {
                std::stringstream ss;
                ss << (size_t)this;
                m_client_id = ss.str();
            }
			
			// Initialise ORB
			int         i            = 0;
			char**      c            = 0;
			const char* options[][2] = { { (char*)"traceLevel", debug_level.c_str()}, { 0, 0 } };
			
			m_orb = CORBA::ORB_init(i, c, "omniORB4", options);
			
			// Bind ORB object to name service object
			CORBA::Object_var obj = m_orb->resolve_initial_references("NameService");
			assert (!CORBA::is_nil(obj.in()));
			
			// Narrow this to the naming context
			CosNaming::NamingContext_var context  = CosNaming::NamingContext::_narrow(obj.in());
			assert (!CORBA::is_nil(context.in()));
			
			// Bind to the name server and lookup 
			CosNaming::Name m_name;
			m_name.length(1);
			m_name[0].id = CORBA::string_dup(service_id.c_str());
			
			// Resolve object reference.
			CORBA::Object_var orid = context->resolve(m_name);
			assert(!CORBA::is_nil(orid.in()));
			
			m_rrsi       = RRSInterface::_narrow(orid.in());
			if (CORBA::is_nil(m_rrsi.in()))
				std::cerr << "IOR is not an SA object reference." << std::endl;
			
		} catch (CORBA::COMM_FAILURE&)        {
			
			std::cerr << "Caught system exception COMM_FAILURE -- unable to contact the object.\n" << std::endl;
			throw DS_ServerConnectionException();
			return;
			
		} catch (CORBA::SystemException&)     {
            
			std::cerr << "Caught a CORBA::SystemException." << std::endl;
			throw DS_SystemException();
			return;
			
		} catch (CORBA::Exception&)           {
			
			std::cerr << "Caught CORBA::Exception." << std::endl;
			throw DS_Exception();
			return;
			
		} catch (omniORB::fatalException& fe) {
			
			std::cerr << "Caught omniORB::fatalException:" << std::endl;
			std::cerr << "  file: " << fe.file() << std::endl;
			std::cerr << "  line: " << fe.line() << std::endl;
			std::cerr << "  mesg: " << fe.errmsg() << std::endl;
			throw DS_FatalException();
			return;
			
		} catch(...) {
			
			std::cerr << "Caught unknown exception." << std::endl;
			throw DS_Exception();
			return;
			
		}
		
	}
	
	
	
	RemoteConnector::~RemoteConnector         ()            {
		m_rrsi->CleanUp();
		m_orb->destroy();
	}
	
	
	
	codeare::error_code
	RemoteConnector::Init (const std::string& name) {
		
		// Prepare configuration for the journey
		std::stringstream  temp;
		temp << GetConfig();
		m_rrsi->config  (temp.str().c_str());
		
		// Initialise back end
        if (m_rrsi->Init (name.c_str(), m_client_id.c_str()) != codeare::OK)
	  return codeare::CANNOT_LOAD_LIBRARY;
		
		return (codeare::error_code)codeare::OK;
		
	}
	
	
	
	codeare::error_code 
	RemoteConnector::Prepare  (const std::string& name)  {
        std::string call = m_client_id + name;
		return (codeare::error_code) m_rrsi->Prepare (call.c_str());
	}
	
	
	
	codeare::error_code 
	RemoteConnector::Process  (const std::string& name)  {
        std::string call = m_client_id + name;
		return  (codeare::error_code) m_rrsi->Process (call.c_str());
	}
	
	
	
	codeare::error_code
	RemoteConnector::Finalise (const std::string& name) {
        std::string call = m_client_id + name;
		return (codeare::error_code) m_rrsi->Finalise(call.c_str());
	}

	
	
	inline long
	RemoteConnector::GetSize (const longs& dims) const {
		
		long size = 1;

		for (int i = 0; i < dims.length(); i++)
			size *= dims[i];

		return size;
		
	}

}
