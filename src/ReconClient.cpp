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

#include "ReconClient.hpp"
#include "CorbaExceptions.hpp"
#include "Matrix.hpp"

using namespace RRClient;

#include <assert.h>
#include <complex>


/**************************************************************************************************/
ReconClient::ReconClient          (const char* name, const char* debug) {

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


/**************************************************************************************************/
ReconClient::~ReconClient         ()            {

	m_orb->destroy();
    
}


/**************************************************************************************************/
error_code
ReconClient::Init (const char* name) {

    error_code  result  = OK;

	// Prepare configuration for the journey
	std::stringstream  temp;
	temp << GetConfig();
	m_rrsi->config  (temp.str().c_str());

	// Initialise back end
    m_rstrats.push_back (m_rrsi->Init (name));
    
	if (m_rstrats.back() == -1)
		result = CANNOT_LOAD_LIBRARY;

	return result;

}


/**************************************************************************************************/
error_code 
ReconClient::Process  (const char* name)  {
    
    error_code  result  = OK;

	// Actually process data
    result    = m_rrsi->Process (m_rstrats.back());
	SetConfig (m_rrsi->config());

    return result;
    
}


/**************************************************************************************************/
error_code
ReconClient::Finalise (const char* name) {

	return m_rrsi->Finalise(m_rstrats.back());

}


/**************************************************************************************************/
long 
ReconClient::GetSize (longs dims) {

    long size = 1;

	for (int i = 0; i < INVALID_DIM; i++)
		size *= dims[i];

    return size;
    
}


/**************************************************************************************************/
void 
ReconClient::SetRaw (Matrix< std::complex<float> >& M) {
	
	raw_data r; 
	
	r.dims.length(INVALID_DIM);
	
	for (int j = 0; j < INVALID_DIM; j++)
		r.dims[j] = M.Dim(j);
	
	r.dreal.length(M.Size()); 
	r.dimag.length(M.Size());
	
	for (int i = 0; i < M.Size(); i++) {
		r.dreal[i] = M[i].real();
		r.dimag[i] = M[i].imag(); 
	}
	
	m_rrsi->raw(r);
	
}
		    

/**************************************************************************************************/
void 
ReconClient::GetRaw (Matrix< std::complex<float> >& M) {
	
	raw_data* rp = m_rrsi->raw();
	
	for (int j = 0; j < INVALID_DIM; j++)
		M.Dim(j) = rp->dims[j];
	
	M.Reset();
	
	for (int i = 0; i < GetSize(rp->dims); i++)
		M[i] = std::complex<float>(rp->dreal[i],rp->dimag[i]);
	
}


/**************************************************************************************************/
void 
ReconClient::SetRHelper (Matrix< std::complex<float> >& M) {
	
	raw_data r; 
	
	r.dims.length(INVALID_DIM);
	
	for (int j = 0; j < INVALID_DIM; j++)
		r.dims[j] = M.Dim(j);
	
	r.dreal.length(M.Size()); 
	r.dimag.length(M.Size());
	
	for (int i = 0; i < M.Size(); i++) {
		r.dreal[i] = M[i].real();
		r.dimag[i] = M[i].imag(); 
	}
	
	m_rrsi->rhelper(r);
	
}
		    

/**************************************************************************************************/
void 
ReconClient::GetRHelper (Matrix< std::complex<float> >& M) {
	
	raw_data* rp = m_rrsi->rhelper();
	
	for (int j = 0; j < INVALID_DIM; j++)
		M.Dim(j) = rp->dims[j];
	
	M.Reset();
	
	for (int i = 0; i < GetSize(rp->dims); i++)
		M[i] = std::complex<float>(rp->dreal[i],rp->dimag[i]);
	
}



void
ReconClient::SetHelper              (Matrix< double >& M) {
	
	helper_data h;
	
	h.dims.length(INVALID_DIM);
	
	for (int j = 0; j < INVALID_DIM; j++)
		h.dims[j] = M.Dim(j);
	
	h.vals.length(M.Size());
	
	for (int i = 0; i < M.Size(); i++)
		h.vals[i] = M[i];
	
	m_rrsi->helper(h);
	
}
		

void
ReconClient::GetHelper              (Matrix< double >& M) {
	
	helper_data* hp = m_rrsi->helper();
	
	for (int j = 0; j < INVALID_DIM; j++)
		M.Dim(j) = hp->dims[j];
	
	M.Reset();
	
	for (int i = 0; i < GetSize(hp->dims); i++)
		M[i] = hp->vals[i];
	
}


void
ReconClient::SetKSpace              (Matrix< double >& M) {
	
	helper_data h;
	
	h.dims.length(INVALID_DIM);
	
	for (int j = 0; j < INVALID_DIM; j++)
		h.dims[j] = M.Dim(j);
	
	h.vals.length(M.Size());
	
	for (int i = 0; i < M.Size(); i++)
		h.vals[i] = M[i];
	
	m_rrsi->kspace(h);
	
}


void
ReconClient::GetKSpace              (Matrix< double >& M) {
			
	helper_data* hp = m_rrsi->kspace();

	for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = hp->dims[j];
	
	M.Reset();
	
	for (int i = 0; i < GetSize(hp->dims); i++)
		M[i] = hp->vals[i];
	
}


void
ReconClient::SetPixel            (Matrix<short>& M) {
			
	pixel_data p;
	
	p.dims.length(INVALID_DIM);
	
	for (int j = 0; j < INVALID_DIM; j++)
		p.dims[j] = M.Dim(j);
	
	p.vals.length(M.Size()); 
	
	for (int i = 0; i < M.Size(); i++)
		p.vals[i] = M[i];
	
	m_rrsi->pixel(p);
	
}


void
ReconClient::GetPixel            (Matrix<short>& M) {
	
	pixel_data* pp = m_rrsi->pixel();

	for (int j = 0; j < INVALID_DIM; j++)
		M.Dim(j) = pp->dims[j];
	
	M.Reset();
	
	for (int i = 0; i < GetSize(pp->dims); i++)
		M[i] = pp->vals[i];
	
}


