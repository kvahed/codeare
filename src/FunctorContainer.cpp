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

#include "FunctorContainer.hpp"

FunctorContainer::~FunctorContainer ()               {
	
	this->Finalise();
	
}


error_code
FunctorContainer::CleanUp () {

	
	this->Finalise();
	return OK;
	
}


error_code
FunctorContainer::Init (const char* name) {
	
	ReconContext* rc;
	
	m_contexts.insert (pair< std::string, ReconContext* > (std::string(name), rc = new ReconContext(name)));

    rc->SetConfig (m_config);

	if ((rc->Init()) != OK) {
		this->Finalise();
		return CONTEXT_CONFIGURATION_FAILED;
	}

	return OK;
	
}


error_code
FunctorContainer::Finalise (const char* name) {
	
	if (!name) {
		
		DataBase::Instance()->Finalise();
		
	} else {
		
		map<string, ReconContext*>::iterator it = m_contexts.find (name);
		
		if (it == m_contexts.end()) {
			Finalise ();
			return CONTEXT_NOT_FOUND;
		} else {
			delete it->second;
			m_contexts.erase(it);
		}
	}
	
	return OK;
	
}


error_code
FunctorContainer::Process  (const char* name)       {
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		return CONTEXT_NOT_FOUND;

	return it->second->Process();
	
}


error_code
FunctorContainer::Prepare  (const char* name)       {
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		return CONTEXT_NOT_FOUND;

	return it->second->Prepare();
	
}



void
FunctorContainer::config (const char* c)    {
		
	std::stringstream tmp;
	
	tmp << c;
	m_config = new char[tmp.str().length() + 1];
	strcpy (m_config, tmp.str().c_str());
		
}
	
