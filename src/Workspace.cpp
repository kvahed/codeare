#include "Workspace.hpp"
#include "BLACS.hpp"

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

		if      (nit->second[1].compare(typeid(cxfl).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<cxfl  > > >(dit->second);
		else if (nit->second[1].compare(typeid(cxdb).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<cxdb  > > >(dit->second);
		else if (nit->second[1].compare(typeid(float).name())  == 0)
			delete boost::any_cast<Ptr<Matrix<float > > >(dit->second);
		else if (nit->second[1].compare(typeid(double).name()) == 0)
			delete boost::any_cast<Ptr<Matrix<double> > >(dit->second);
		else if (nit->second[1].compare(typeid(short).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<short > > >(dit->second);
		else if (nit->second[1].compare(typeid(long).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<long  > > >(dit->second);

  		m_store.erase(dit);
  		m_ref.erase(nit);

  	}

	return OK;
	
}

