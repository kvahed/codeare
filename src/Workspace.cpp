#include "Workspace.hpp"
#include "Algos.hpp"
#include "Print.hpp"


Workspace* Workspace::m_inst = 0; 


Workspace::Workspace() {}


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

codeare::error_code
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
    
	return codeare::OK;
	
}

void Workspace::Print (std::ostream& os) const {

	os << "codeare service " << VERSION << endl;
#ifdef GIT_COMMIT
	os << "commit " << GIT_COMMIT << " [" << GIT_COMMIT_DATE << "]" << endl;
#endif

	os << "\nMatrices\n";
	os <<   "--------\n\n";
	for (reflist::const_iterator i = m_ref.begin(); i != m_ref.end(); i++) {

	    const boost::any& b = m_store.find (i->second[0])->second;
	    const std::string k_name = i->first;
	    const size_t kl = k_name.length();

	    os << setw(24) << k_name << " | ";
	    if (b.type() == typeid(boost::shared_ptr<Matrix<float> >))
		    os << "           float |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<float> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<double> >))
		    os << "          double |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<double> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<cxfl> >))
		    os << "  complex<float> |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<cxfl> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<cxdb> >))
		    os << " complex<double> |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<cxdb> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<short> >))
		    os << "           short |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<short> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<long> >))
		    os << "            long |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<long> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<size_t> >))
		    os << "          size_t |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<size_t> > >(b));
	    else if (b.type() == typeid(boost::shared_ptr<Matrix<cbool> >))
		    os << "            bool |" << setw(8) << size(*boost::any_cast<boost::shared_ptr<Matrix<cbool> > >(b));
	}

    os << "\n\nParameters\n" ;
    os <<       "----------\n\n";

    os << p;

}



