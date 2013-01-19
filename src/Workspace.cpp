#include "Workspace.hpp"


Workspace* Workspace::m_inst = 0; 

Workspace::~Workspace () { 
	
	Finalise();
	m_inst = 0; 
	
}


Workspace* 
Workspace::Instance ()  {

    if (m_inst == 0) 
        m_inst = new Workspace ();

	bool t = true;
	bool f = false;
	
	m_inst->SetAttribute("FFTWFThreadsInitialised", &f);
	m_inst->SetAttribute("FFTWThreadsInitialised", &f);
	
    return m_inst;
		
}
	
error_code
Workspace::Finalise () {
	
  	while (!m_store.empty())
		m_store.erase(m_store.begin());
			   
  	while (!m_ref.empty())
		m_ref.erase(m_ref.begin());

	return OK;
	
}

