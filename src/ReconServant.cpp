/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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
#include "ReconContext.hpp"

using namespace RRServer;

ReconServant::ReconServant  () : 
	m_config (new char) {

}


ReconServant::~ReconServant ()               {

	this->CleanUp();

}


error_code
ReconServant::CleanUp () {

	this->Finalise();
	delete m_config;

}

error_code
ReconServant::Init (const char* name) {

	ReconContext* rc;
	error_code e = OK;
	
	m_contexts.insert (pair< std::string, ReconContext* > (std::string(name), rc = new ReconContext(name)));

    rc->SetConfig (m_config);

	if ((e = rc->Init()) != OK) {
		Finalise();
		return CONTEXT_CONFIGURATION_FAILED;
	}

	return e;

}


error_code
ReconServant::Finalise (const char* name) {

	error_code e = OK;

	if (!name) {

		DataBase::Instance()->Finalise();

	} else {
		
		map<string, ReconContext*>::iterator it = m_contexts.find (name);
		
		if (it == m_contexts.end()) {
			e = CONTEXT_NOT_FOUND;
		} else {
			delete it->second;
			m_contexts.erase(it);
		}
	}

	return e;

}
	

error_code
ReconServant::Process  (const char* name)       {
	
	error_code e = OK;
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		e = CONTEXT_NOT_FOUND;
	else {
		std::cout << "Processing ... " << name;
		error_code e = it->second->Process();
		std::cout << " ... done. " << std::endl;
	}
	
	return e;
	
}


error_code
ReconServant::Prepare  (const char* name)       {
	
	error_code e = OK;
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		e = CONTEXT_NOT_FOUND;
	else {
		std::cout << "Preparing ... " << std::endl;
		error_code e = it->second->Prepare();
		std::cout << "... done. " << std::endl;
	}
	
	return e;
	
}


void
ReconServant::set_cxfl  (const char* name, const cxfl_data& c) {
	
	DataBase::Instance()->SetCXFL(name, c);
	
}


void
ReconServant::get_cxfl (const char* name, cxfl_data& c) {
	
	c.dims.length(INVALID_DIM);
	c.res.length (INVALID_DIM);
	DataBase::Instance()->GetCXFL(name, c);
	
}


void
ReconServant::set_cxdb  (const char* name, const cxdb_data& c) {
	
	DataBase::Instance()->SetCXDB(name, c);
	
}


void
ReconServant::get_cxdb (const char* name, cxdb_data& c) {
	
	c.dims.length(INVALID_DIM);
	c.res.length (INVALID_DIM);
	DataBase::Instance()->GetCXDB(name, c);
	
}


void
ReconServant::set_rldb       (const char* name, const rldb_data& r)   {
	
	DataBase::Instance()->SetRLDB(name, r);
	
}


void
ReconServant::get_rldb       (const char* name, rldb_data& r) {
	
	r.dims.length(INVALID_DIM);
	r.res.length (INVALID_DIM);
	DataBase::Instance()->GetRLDB(name, r);
	
}


void
ReconServant::set_rlfl       (const char* name, const rlfl_data& r)   {
	
	DataBase::Instance()->SetRLFL(name, r);
	
}


void
ReconServant::get_rlfl       (const char* name, rlfl_data& r) {
	
	r.dims.length(INVALID_DIM);
	r.res.length(INVALID_DIM);
	DataBase::Instance()->GetRLFL(name, r);

}


void
ReconServant::set_shrt        (const char* name, const shrt_data& p) {

	DataBase::Instance()->SetSHRT(name, p);

}


void
ReconServant::get_shrt        (const char* name, shrt_data& p) {

	p.dims.length(INVALID_DIM);
	p.res.length(INVALID_DIM);
	DataBase::Instance()->GetSHRT(name, p);

}


void
ReconServant::set_long        (const char* name, const long_data& p) {

	DataBase::Instance()->SetLONG(name, p);

}


void
ReconServant::get_long        (const char* name, long_data& p) {

	p.dims.length(INVALID_DIM);
	p.res.length(INVALID_DIM);
	DataBase::Instance()->GetLONG(name, p);

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

