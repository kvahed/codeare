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

#include "RemoteConnector.hpp"
#include "CorbaExceptions.hpp"


#include <assert.h>
#include <complex>

namespace RRClient {


	RemoteConnector::RemoteConnector          (const char* name, const char* debug) {
		
		try {
			
			// Initialise ORB
			int         i            = 0;
			char**      c            = 0;
			const char* options[][2] = { { (char*)"traceLevel", debug }, { 0, 0 } };
			
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
			m_name[0].id = CORBA::string_dup(name);
			
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
	
	
	
	error_code
	RemoteConnector::Init (const char* name) {
		
		// Prepare configuration for the journey
		std::stringstream  temp;
		temp << GetConfig();
		m_rrsi->config  (temp.str().c_str());
		
		// Initialise back end
		if (m_rrsi->Init (name) != OK)
			return CANNOT_LOAD_LIBRARY;
		
		return OK;
		
	}
	
	
	
	error_code 
	RemoteConnector::Prepare  (const char* name)  {
		
		return m_rrsi->Prepare (name);
		
	}
	
	
	
	error_code 
	RemoteConnector::Process  (const char* name)  {
		
		return  m_rrsi->Process (name);
		
	}
	
	
	
	error_code
	RemoteConnector::Finalise (const char* name) {
		
		return m_rrsi->Finalise(name);
		
	}

	
	
	long 
	RemoteConnector::GetSize (const longs dims) const {
		
		long size = 1;
		
		for (int i = 0; i < INVALID_DIM; i++)
			size *= dims[i];
		
		return size;
		
	}
	
	
	
	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix< std::complex<float> >& M) const {
		
		cxfl_data r; 
		
		r.dims.length(INVALID_DIM);
		r.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) {
			r.dims[j] = M.Dim(j);
			r.res[j]  = M.Res(j);
		}
		
		r.vals.length(2 * M.Size());
		
		memcpy (&r.vals[0], &M[0], r.vals.length() * sizeof(float));
		
		m_rrsi->set_cxfl(name.c_str(), r);
		
	}
	
	
	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix< std::complex<float> >& m) const {
		
		cxfl_data rp; m_rrsi->get_cxfl(name.c_str(), rp);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = rp.dims[j];
			mress[j] = rp.res[j]; 
		}
		
		m = Matrix<cxfl> (mdims, mress);
		
		memcpy (&m[0], &rp.vals[0], rp.vals.length() * sizeof(float));
		
	}
	

	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix<cxdb>& M) const {
		
		cxdb_data r; 
		
		r.dims.length(INVALID_DIM);
		r.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) {
			r.dims[j] = M.Dim(j);
			r.res[j]  = M.Res(j);
		}
		
		r.vals.length(2 * M.Size());
		
		memcpy (&r.vals[0], &M[0], r.vals.length() * sizeof(double));
		
		m_rrsi->set_cxdb(name.c_str(), r);
		
	}
	
	
	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix<cxdb>& m) const {
		
		cxdb_data rp; m_rrsi->get_cxdb(name.c_str(), rp);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = rp.dims[j];
			mress[j] = rp.res[j]; 
		}
		
		m = Matrix<cxdb> (mdims, mress);
		
		memcpy (&m[0], &rp.vals[0], rp.vals.length() * sizeof(double));
		
	}
	

	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix<double>& m) const {
		
		rldb_data r;
		
		r.dims.length(INVALID_DIM); r.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) { r.dims[j] = m.Dim(j);
			r.res[j] = m.Res(j); }
		
		r.vals.length(m.Size());
		
		for (int i = 0; i < m.Size(); i++) r.vals[i] = m[i];
		
		m_rrsi->set_rldb(name.c_str(), r);
		
	}
	
	
	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix<double>& m) const {
		
		rldb_data r; m_rrsi->get_rldb(name.c_str(), r);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = r.dims[j];
			mress[j] = r.res[j]; 
		}
		
		m = Matrix<double> (mdims, mress);
		
		for (int i = 0; i < GetSize(r.dims); i++) m[i] = r.vals[i];
		
	}
	
	
	
	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix<float>& m) const {
		
		rlfl_data r;
		
		r.dims.length(INVALID_DIM); r.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) { r.dims[j] = m.Dim(j);
			r.res[j] = m.Res(j); }
		
		r.vals.length(m.Size());
		
		for (int i = 0; i < m.Size(); i++) r.vals[i] = m[i];
		
		m_rrsi->set_rlfl(name.c_str(), r);
		
	}
	

	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix<float>& m) const {
		
		rlfl_data r; m_rrsi->get_rlfl(name.c_str(), r);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = r.dims[j];
			mress[j] = r.res[j]; 
		}
		
		m = Matrix<float> (mdims, mress);
		
		for (int i = 0; i < GetSize(r.dims); i++) m[i] = r.vals[i];
		
	}
	
	
	
	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix<short>& m) const {
		
		shrt_data p;
		
		p.dims.length(INVALID_DIM); p.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) { p.dims[j] = m.Dim(j);
			p.res[j] = m.Res(j); }
		
		p.vals.length(m.Size());
		
		for (int i = 0; i < m.Size(); i++) p.vals[i] = m[i];
		
		m_rrsi->set_shrt(name.c_str(), p);
		
	}
	
	
	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix<short>& m) const {
		
		shrt_data p; m_rrsi->get_shrt(name.c_str(), p);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = p.dims[j];
			mress[j] = p.res[j]; 
		}
		
		m = Matrix<short> (mdims, mress);
		
		for (int i = 0; i < GetSize(p.dims); i++) m[i] = p.vals[i];
		
	}


	template<> void 
	RemoteConnector::SetMatrix (const std::string& name, Matrix<long>& m) const {
		
		long_data p;
		
		p.dims.length(INVALID_DIM); p.res.length(INVALID_DIM);
		
		for (int j = 0; j < INVALID_DIM; j++) { p.dims[j] = m.Dim(j);
			p.res[j] = m.Res(j); }
		
		p.vals.length(m.Size());
		
		for (int i = 0; i < m.Size(); i++) p.vals[i] = m[i];
		
		m_rrsi->set_long(name.c_str(), p);
		
	}
	
	
	template<> void 
	RemoteConnector::GetMatrix (const std::string& name, Matrix<long>& m) const {
		
		long_data p; m_rrsi->get_long(name.c_str(), p);
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];
		
		for (int j = 0; j < INVALID_DIM; j++) { 
			mdims[j] = p.dims[j];
			mress[j] = p.res[j]; 
		}
		
		m = Matrix<long> (mdims, mress);		

		for (int i = 0; i < GetSize(p.dims); i++) m[i] = p.vals[i];
		
	}


}
