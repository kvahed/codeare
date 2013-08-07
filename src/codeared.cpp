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

#include "codeared.h"
#include "MongooseService.hpp"

int main (int argc, char** argv) {

    if (!init (argc, argv))
        return 1;

    // --------------------------------------------------------------------------
    // Start HTTP server:
    // --------------------------------------------------------------------------
    using namespace codeare::service;
    MongooseService& mg = MongooseService::Instance();

    // --------------------------------------------------------------------------
    // Start CORBA server:
    // --------------------------------------------------------------------------
    try {
        
        streambuf* out;
        streambuf* err;
        ofstream   log (logfile);
        
        if (log.is_open()) {
            out = cout.rdbuf(log.rdbuf());
            err = cerr.rdbuf(log.rdbuf());
        } else {
            cout << "Could not open logfile " << logfile << "." << endl;
            cout << "Exiting :(" << endl << endl;
            return 1;
        }
        
        // Initialise ORB
        const char*    options[][2] = { { (char*)"traceLevel", debug}, /*{ (char*)"traceFile", logfile}, */{ 0, 0 } };
        CORBA::ORB_var orb          = CORBA::ORB_init(argc, argv, "omniORB4", options);
        
        // Get reference to the RootPOA.
        CORBA::Object_var obj = orb->resolve_initial_references("RootPOA");
        PortableServer::POA_var _poa = PortableServer::POA::_narrow(obj.in());
        
        // Initialise servant
        ReconServant* myReconServant = new ReconServant();

        // Activate in RootPDA
        PortableServer::ObjectId_var myReconServant_oid
            = _poa->activate_object(myReconServant);
        
        // Obtain object reference from servant and register in naming service
        CORBA::Object_var SA_obj = myReconServant->_this();
        
        // Obtain a reference to the object, and print it out as string IOR.
        CORBA::String_var sior(orb->object_to_string(SA_obj.in()));
        
        // Bind to the name server and lookup 
        CORBA::Object_var obj1=orb->resolve_initial_references("NameService");
        assert(!CORBA::is_nil(obj1.in()));
        
        // Get context
        CosNaming::NamingContext_var nc = CosNaming::NamingContext::_narrow(obj1.in());
        assert(!CORBA::is_nil(nc.in()));
        
        // Resolve name
        CosNaming::Name m_name;
        m_name.length(1);
        m_name[0].id=CORBA::string_dup(name);
        nc->rebind (m_name,SA_obj.in());
        
        // Activate POA manager
        PortableServer::POAManager_var pmgr = _poa->the_POAManager();
        pmgr->activate();
        
        // Accept requests from clients
        orb->run();
                                                                                
        orb->destroy();
        
        free(m_name[0].id); // str_dup does a malloc internally
        
        cout.rdbuf(out);
        cerr.rdbuf(err);
        
    } catch(CORBA::SystemException&) {
        cerr << "Caught CORBA::SystemException." << endl;
        throw DS_SystemException();
    } catch(CORBA::Exception&) {
        cerr << "Caught CORBA::Exception." << endl;
        throw DS_Exception();
    } catch(omniORB::fatalException& fe) {
        cerr << "Caught omniORB::fatalException:" << endl;
        cerr << "  file: " << fe.file() << endl;
        cerr << "  line: " << fe.line() << endl;
        cerr << "  mesg: " << fe.errmsg() << endl;
        throw DS_FatalException();
    } catch(...) {
        cerr << "Caught unknown exception." << endl;
        throw DS_Exception();
    }
    
    return 0;
    
}


