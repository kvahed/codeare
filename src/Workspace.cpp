#include "Workspace.hpp"

Workspace* Workspace::m_inst = 0; 


Workspace::Workspace() {
    
	bool f = false;
	int  zero = 0;
    
	p["FFTWFThreads"] = zero;
	p[ "FFTWThreads"] = zero;
    
}


Workspace::~Workspace () { 
    
	Finalise();
	m_inst = 0;
    
}


Workspace&
Workspace::Instance ()  {

    if (m_inst == 0)
        m_inst = new Workspace ();

	return *m_inst;
	
}

error_code
Workspace::Finalise () {
    
	while (!m_ref.empty()) {
        
  		reflist::iterator nit = m_ref.begin();
  		store::iterator dit = m_store.find (nit->second[0]);
        
		if      (nit->second[1].compare(typeid(cxfl).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<cxfl > > >(dit->second).reset();
		else if (nit->second[1].compare(typeid(cxdb).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<cxdb > > >(dit->second).reset();
		else if (nit->second[1].compare(typeid(float).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<float > > >(dit->second).reset();
		else if (nit->second[1].compare(typeid(double).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<double> > >(dit->second).reset();
		else if (nit->second[1].compare(typeid(short).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<short > > >(dit->second).reset();
		else if (nit->second[1].compare(typeid(long).name()) == 0)
			boost::any_cast<boost::shared_ptr<Matrix<long > > >(dit->second).reset();
        
  		m_store.erase(dit);
  		m_ref.erase(nit);
        
  	}
    
	return OK;
	
}

