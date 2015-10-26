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

Queue::~Queue () {}


short Queue::CleanUp () {
	this->Finalise();
	return (short) codeare::OK;
}


short Queue::Init (const char* name, const char* config, const char* client_id) {
	ReconContext* rc = new ReconContext(name);
    rc->SetConfig (config);
	if ((rc->Init()) != codeare::OK) {
		this->Finalise();
		return (short) codeare::CONTEXT_CONFIGURATION_FAILED;
	}
	m_contexts.push_back (QEntry(std::string(client_id) + std::string(name), rc));
	return (short) codeare::OK;
}


short Queue::Finalise (const char* name) {
	while (!m_contexts.empty()) {
		auto it = m_contexts.begin();
		delete it->context;
		m_contexts.erase(it);
	}
	Workspace::Instance().Finalise();
	return (short) codeare::OK;
}


short Queue::Process  (const char* name)       {
    codeare::error_code ret = codeare::OK;
	for (auto it = m_contexts.begin(); it != m_contexts.end(); ++it) {
		cout << it->name << endl;
		SimpleTimer t(name);
		ret = it->context->Process();
		t.Stop();
		if (ret != codeare::OK) {
			printf ("Procession of %s failed\n", it->name.c_str());
			break;
		}
	}
	return (short)ret;
}


short Queue::Prepare  (const char* name)       {
	short ret = 0;
	for (auto it = m_contexts.begin(); it != m_contexts.end(); ++it)
		if ((ret = it->context->Prepare()) != codeare::OK) {
			printf ("Preparation of %s \n", it->name.c_str());
			break;
		}
	return (short)ret;
}


void Queue::config (const char* c)    {
	std::stringstream tmp;
	tmp << c;
	m_config = new char[tmp.str().length() + 1];
	strcpy (m_config, tmp.str().c_str());
}
	
