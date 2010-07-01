#include "ReconClient.h"

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
		cerr << "IOR is not an SA object reference." << endl;

	} catch (CORBA::COMM_FAILURE&)        {

		cerr << "Caught system exception COMM_FAILURE -- unable to contact the object.\n" << endl;
		throw DS_ServerConnectionException();
		return;
		
	} catch (CORBA::SystemException&)     {
		
		cerr << "Caught a CORBA::SystemException." << endl;
		throw DS_SystemException();
		return;
		
	} catch (CORBA::Exception&)           {
		
		cerr << "Caught CORBA::Exception." << endl;
		throw DS_Exception();
		return;
		
	} catch (omniORB::fatalException& fe) {

		cerr << "Caught omniORB::fatalException:" << endl;
		cerr << "  file: " << fe.file() << endl;
		cerr << "  line: " << fe.line() << endl;
		cerr << "  mesg: " << fe.errmsg() << endl;
		throw DS_FatalException();
		return;
		
	} catch(...) {
		
		cerr << "Caught unknown exception." << endl;
		throw DS_Exception();
		return;
		
	}
	
	m_raw    = new raw_data;
	m_helper = new raw_data;
	m_pixel  = new pixel_data;
	
	m_raw->dims.length(INVALID_DIM);
	m_helper->dims.length(INVALID_DIM);
	m_pixel->dims.length(INVALID_DIM);
	
	return;
	
};



/**************************************************************************************************/
ReconClient::~ReconClient         ()            {

	delete m_raw;
	delete m_helper;
	delete m_pixel;

	m_orb->destroy();
	
};



/**************************************************************************************************/
error_code 
ReconClient::Requestprocess_data  (method m)  {
	
	error_code result = OK;
	
	// Set data for recon
	m_rrsi->raw(m_raw[0]);
	
	// Reconstruct on remote service
	result = m_rrsi->process_data(m);
	
	// Get data back from recon
	if (result == OK) {
		
		if(m_have_raw == true)
			m_raw = m_rrsi->raw();
		if (m_have_helper == true)
			m_helper = m_rrsi->helper();
		if (m_have_pixel == true)
			m_pixel = m_rrsi->pixel();
		
		m_labels = m_rrsi->labels();
		
	}
	
	return result;
	
};



/**************************************************************************************************/
long 
ReconClient::GetRawSize () {

	long size = 1;
	for (int i = 0; i < INVALID_DIM; i++)
		size *= m_raw->dims[i];
	return size;
	
};



/**************************************************************************************************/
long 
ReconClient::GetHelperSize () {
	
	long size = 1;
	for (int i = 0; i < INVALID_DIM; i++)
		size *= m_helper->dims[i];
	return size;
	
};



/**************************************************************************************************/
long 
ReconClient::GetPixelSize () {
	
	long size = 1;
	for (int i = 0; i < INVALID_DIM; i++)
		size *= m_pixel->dims[i];
	return size;
	
};

	
