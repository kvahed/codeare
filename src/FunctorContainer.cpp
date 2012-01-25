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
	error_code e = OK;
	
	m_contexts.insert (pair< std::string, ReconContext* > (std::string(name), rc = new ReconContext(name)));
	
    rc->SetConfig (m_config);
	
	if ((e = rc->Init()) != OK) {
		this->Finalise();
		return CONTEXT_CONFIGURATION_FAILED;
	}
	
	return e;
	
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
FunctorContainer::Prepare  (const char* name)       {
	
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



