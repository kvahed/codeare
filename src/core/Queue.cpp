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

#include "Queue.hpp"

Queue::~Queue ()               {
	
	this->Finalise();
	
}


short Queue::CleanUp () {

	this->Finalise();
	return (short) codeare::OK;
	
}


short Queue::Init (const char* name, const char* client_id) {
	
	ReconContext* rc;
	
	m_contexts.insert (pair< std::string, ReconContext* > (std::string(client_id) + std::string(name), rc = new ReconContext(name)));

    rc->SetConfig (m_config);

	if ((rc->Init()) != codeare::OK) {
		this->Finalise();
		return codeare::CONTEXT_CONFIGURATION_FAILED;
	}

	return (short) codeare::OK;
	
}


short Queue::Finalise (const char* name) {
	
	if (!name) {
		
		Workspace::Instance().Finalise();
		
	} else {
		
		map<string, ReconContext*>::iterator it = m_contexts.find (name);
		
		if (it == m_contexts.end()) {
			Finalise ();
			return codeare::CONTEXT_NOT_FOUND;
		} else {
			delete it->second;
			m_contexts.erase(it);
		}
	}
	
	return (short) codeare::OK;
	
}


short Queue::Process  (const char* name)       {
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		return codeare::CONTEXT_NOT_FOUND;

	return (short) it->second->Process();
	
}


short Queue::Prepare  (const char* name)       {
	
	map<string, ReconContext*>::iterator it = m_contexts.find (name);
	
	if (it == m_contexts.end()) 
		return codeare::CONTEXT_NOT_FOUND;

	return (short) it->second->Prepare();
	
}



void Queue::config (const char* c)    {
		
	std::stringstream tmp;
	
	tmp << c;
	m_config = new char[tmp.str().length() + 1];
	strcpy (m_config, tmp.str().c_str());
		
}
	
