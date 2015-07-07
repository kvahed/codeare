#include "Workspace.hpp"
#include "Algos.hpp"
#include "Print.hpp"


Workspace* Workspace::m_inst = 0; 

Workspace::Workspace () {}

Workspace::~Workspace () { 
	Finalise();
	m_inst = 0;
}


codeare::error_code
Workspace::Finalise () {
    

	while (!m_ref.empty()) {
  		reflist::iterator nit = m_ref.begin();
  		m_store.erase(m_store.find (nit->second[0]));
  		m_ref.erase(nit);
  	}
    
	return codeare::OK;
	
}

void Workspace::Print (std::ostream& os) const {

    os << "codeare service ";
#ifdef VERSION
	os << VERSION;
#endif
    os << endl;
#ifdef GIT_COMMIT
	os << "commit " << GIT_COMMIT << " [" << GIT_COMMIT_DATE << "]" << endl;
#endif

	os << "\nMatrices\n";
	os <<   "--------\n\n";
	for (reflist::const_iterator i = m_ref.begin(); i != m_ref.end(); i++) {

	    const boost::any& b = m_store.find (i->second[0])->second;
	    const std::string k_name = i->first;
	    const size_t kl = k_name.length();
#ifndef __INTEL_COMPILER
	    os << setw(24) << k_name << " | ";
	    if (b.type() == typeid(shrd_ptr<Matrix<float> >))
		    os << "           float |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<float> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<double> >))
		    os << "          double |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<double> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<cxfl> >))
		    os << "  complex<float> |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<cxfl> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<cxdb> >))
		    os << " complex<double> |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<cxdb> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<short> >))
		    os << "           short |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<short> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<long> >))
		    os << "            long |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<long> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<size_t> >))
		    os << "          size_t |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<size_t> > >(b));
	    else if (b.type() == typeid(shrd_ptr<Matrix<cbool> >))
		    os << "            bool |" << setw(8) << size(*boost::any_cast<shrd_ptr<Matrix<cbool> > >(b));
#endif
	}
    os << "\n\nParameters\n" ;
    os <<       "----------\n\n";

    os << p;

}



