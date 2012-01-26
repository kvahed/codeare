#include "FunctorContainer.hpp"

FunctorContainer::~FunctorContainer ()               {
	
	this->Finalise();
	
}


error_code
FunctorContainer::CleanUp () {

	
	this->Finalise();
	
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
	
