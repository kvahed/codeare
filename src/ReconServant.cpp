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

#include "ReconServant.hpp"
#include "Workspace.hpp"

using namespace RRStrategy;

namespace RRServer {


	ReconServant::ReconServant  () {}
	
	ReconServant::~ReconServant () {}
	
    short
	ReconServant::CleanUp () {
		return Queue::CleanUp();
	}

	void
	ReconServant::config       (const char* d)    {
		std::stringstream tmp;
		tmp << d;
		m_config = new char[tmp.str().length() + 1];
		strcpy (m_config, tmp.str().c_str());
	}
	
	char*
	ReconServant::config       ()                    {
		std::stringstream tmp;
		tmp << m_config;
		return CORBA::string_dup(tmp.str().c_str());
	}
	
	short
	ReconServant::Init (const char* name, const char* config, const char* client_id) {
		return Queue::Init(name, config, client_id);
	}
	
	short
	ReconServant::Finalise (const char* name) {
		return Queue::Finalise(name);
	}
	
    short
	ReconServant::Process  (const char* name) {
		return Queue::Process(name);
	}
	
	short
	ReconServant::Prepare  (const char* name)       {
		return Queue::Prepare(name);
	}
	
	void
	ReconServant::set_cxfl  (const char* name, const RRSModule::cxfl_data& c) {
		this->SetMatrix (name,c);
	}
	
	void
	ReconServant::get_cxfl (const char* name, RRSModule::cxfl_data& c) {
		this->GetMatrix (name, c);
	}

	void
	ReconServant::set_cxdb  (const char* name, const RRSModule::cxdb_data& c) {
		this->SetMatrix (name,c);
	}

	void
	ReconServant::get_cxdb (const char* name, RRSModule::cxdb_data& c) {
		this->GetMatrix (name,c);
	}

	void
	ReconServant::set_rlfl  (const char* name, const RRSModule::rlfl_data& c) {
		this->SetMatrix (name,c);
	}
	
	void
	ReconServant::get_rlfl (const char* name, RRSModule::rlfl_data& c) {
		this->SetMatrix (name,c);
	}
	
	void
	ReconServant::set_rldb  (const char* name, const RRSModule::rldb_data& c) {
		this->SetMatrix (name,c);
	}

	void
	ReconServant::get_rldb (const char* name, RRSModule::rldb_data& c) {
		this->GetMatrix (name,c);
	}

	void
	ReconServant::set_shrt  (const char* name, const RRSModule::shrt_data& c) {
		this->SetMatrix (name,c);
	}

	void
	ReconServant::get_shrt (const char* name, RRSModule::shrt_data& c) {
		this->GetMatrix (name,c);
	}
	
	void
	ReconServant::set_long  (const char* name, const RRSModule::long_data& c) {
		this->SetMatrix (name,c);
	}

	void
	ReconServant::get_long (const char* name, RRSModule::long_data& c) {
		this->GetMatrix (name,c);
	}


    void
    ReconServant::inform (omni::omniInterceptors::assignUpcallThread_T::info_T &info) {
        info.run();
    }

}
