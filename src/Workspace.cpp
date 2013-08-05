#include "Workspace.hpp"
#include "Algos.hpp"
#include "Print.hpp"


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

void Workspace::Print (std::ostream& os) const {

       os << "\nMatrices:\n";
       for (reflist::const_iterator i = m_ref.begin(); i != m_ref.end(); i++) {

           const boost::any& b = m_store.find (i->second[0])->second;
           os << i->first << "\t";
           if (b.type() == typeid(boost::shared_ptr<Matrix<float> >))
        	   os << "\t Matrix<float>          \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<float> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<double> >))
        	   os << "\t Matrix<double>         \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<double> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<cxfl> >))
        	   os << "\t Matrix<complex<float>> \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<cxfl> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<cxdb> >))
        	   os << "\t Matrix<complex<double>>\t" << size(*boost::any_cast<boost::shared_ptr<Matrix<cxdb> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<short> >))
        	   os << "\t Matrix<short>          \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<short> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<long> >))
        	   os << "\t Matrix<long>           \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<long> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<size_t> >))
        	   os << "\t Matrix<size_t>         \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<size_t> > >(b));
           else if (b.type() == typeid(boost::shared_ptr<Matrix<cbool> >))
        	   os << "\t Matrix<bool>           \t" << size(*boost::any_cast<boost::shared_ptr<Matrix<cbool> > >(b));
       }

       os << "\nParameters:\n" ;
       os << p;

   }



