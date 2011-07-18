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
    
    m_raw     = new raw_data;
    m_rhelper = new raw_data;
    m_helper  = new helper_data;
    m_kspace  = new helper_data;
    m_pixel   = new pixel_data;

    m_raw->dims.length(INVALID_DIM);
    m_rhelper->dims.length(INVALID_DIM);
    m_helper->dims.length(INVALID_DIM);
    m_kspace->dims.length(INVALID_DIM);
    m_pixel->dims.length(INVALID_DIM);

	for (int i = 0; i < INVALID_DIM; i++) {
		m_raw->dims[i] = 0;
		m_rhelper->dims[i] = 0;
		m_helper->dims[i] = 0;
		m_kspace->dims[i] = 0;
		m_pixel->dims[i] = 0;
	}

}



/**************************************************************************************************/
ReconClient::~ReconClient         ()            {

    delete m_raw;
    delete m_rhelper;
    delete m_helper;
    delete m_kspace;
    delete m_pixel;

    m_orb->destroy();
    
}


error_code
ReconClient::Init (const char* name) {

    error_code  result  = OK;

    m_rstrats.push_back (m_rrsi->Init (name));
    
	if (m_rstrats.back() == -1)
		result = CANNOT_LOAD_LIBRARY;

	return result;

}


/**************************************************************************************************/
error_code 
ReconClient::Process  (const char* name)  {
    
    error_code  result  = OK;
	std::string confstr = "";

	// Prepare configuration for the journey
	std::stringstream  temp;
	temp << GetConfig();

    // Send data for procession to remote interface
    m_rrsi->raw     (*(m_raw));
    m_rrsi->rhelper (*(m_rhelper));
    m_rrsi->helper  (*(m_helper));
    m_rrsi->kspace  (*(m_kspace));
    m_rrsi->pixel   (*(m_pixel));
	m_rrsi->config  (temp.str().c_str());
	
	delete m_raw;
	delete m_rhelper;
	delete m_helper;
	delete m_kspace;
	delete m_pixel;

    result    = m_rrsi->Process (m_rstrats.back());
    
    m_raw     = new raw_data;
	m_rhelper = new raw_data;
    m_helper  = new helper_data;
    m_kspace  = new helper_data;
    m_pixel   = new pixel_data;

    // Get data back from remote interface
	m_raw     = m_rrsi->raw();

	m_rhelper = m_rrsi->rhelper();
	m_helper  = m_rrsi->helper();
	m_kspace  = m_rrsi->kspace();
	m_pixel   = m_rrsi->pixel();
	SetConfig (m_rrsi->config());

	result    = m_rrsi->Finalise(m_rstrats.back());

    return result;
    
}



/**************************************************************************************************/
long 
ReconClient::GetRawSize () {

    long size = 1;
    for (int i = 0; i < INVALID_DIM; i++)
        size *= m_raw->dims[i];
    return size;
    
}



/**************************************************************************************************/
long 
ReconClient::GetRHelperSize () {
    
    long size = 1;
    for (int i = 0; i < INVALID_DIM; i++)
        size *= m_rhelper->dims[i];
    return size;
    
}



/**************************************************************************************************/
long 
ReconClient::GetHelperSize () {
    
    long size = 1;
    for (int i = 0; i < INVALID_DIM; i++)
        size *= m_helper->dims[i];
    return size;
    
}



/**************************************************************************************************/
long 
ReconClient::GetKSpaceSize () {
    
    long size = 1;
    for (int i = 0; i < INVALID_DIM; i++)
        size *= m_kspace->dims[i];
    return size;
    
}



/**************************************************************************************************/
long 
ReconClient::GetPixelSize () {
    
    long size = 1;
    for (int i = 0; i < INVALID_DIM; i++)
        size *= m_pixel->dims[i];
    return size;
    
}

    
