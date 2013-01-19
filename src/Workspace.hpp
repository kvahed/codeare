#ifndef __WORK_SPACE_HPP_
#define __WORK_SPACE_HPP_

#include "Matrix.hpp"
#include "Configurable.hpp"

#include <boost/any.hpp>
#include <map>

#ifdef __WIN32__ 
  #include "RRSModule.h"
#else
  #include "RRSModule.hh"
#endif

#ifdef __APPLE__
  #include "AppleDigest.hpp"
#else
  #include "Digest.hpp"
#endif

using namespace std;
using namespace RRSModule;


/**
 * @brief Central database for all shared matrices<br/>
 *        Matrix smart pointers are stored in individual (type) maps along with string identifier and retrieved by their names.
 */
class Workspace : public Configurable {


	typedef map<string, string[2]> reflist;
	typedef pair<string, string[2]> ref;
	typedef map<string, boost::any> store;
	typedef pair<string, boost::any> entry;

 public:


	/**
	 * @brief        Clean up
	 */
	virtual ~Workspace        ();


	/**
	 * @brief        Get pointer to database instance
	 */
	static Workspace* 
	Instance         () ;
	

	/**
	 * @brief        Initialise database
	 */
	error_code 
	Initialise       ();


	/**
	 * @brief        Free RAM
	 *
	 * @return       Success
	 */ 
	error_code 
	Finalise         ();
	
	
	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> void
	GetMatrix          (const string name, Matrix<T>& m) {

		if (m_ref.find (name) == m_ref.end())
				return;

		reflist::iterator it = m_ref.find(name);
		m = *boost::any_cast<Ptr<Matrix<T> > >(m_store[it->second[0]]);

	}
	
	
	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> void
	SetMatrix          (const string name, Matrix<T>& m) {

		string tag[2];
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();

		Ptr<Matrix<T> > pm = NEW (Matrix<T>());
		boost::any val     = pm;

		if (m_ref.find (name) == m_ref.end()) {
			m_ref.insert (ref(name,tag));
			m_store.insert (entry(tag[0], val));
		}

		m.SetClassName(name.c_str());

		*pm = m;

	}
	
	
	/**
	 * @brief        Add a matrix to according container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	template <class T> Matrix<T>& 
	AddMatrix        (const string name, Ptr< Matrix<T> > m) {

		std::string tag[2];
		tag[1] = sha256(name);
		tag[0] = typeid(T).name();

		boost::any value = m;

		assert (m_ref.find (name) == m_ref.end());
		m_ref.insert (pair<string, string[2]> (name, tag));
		m_store.insert (entry (tag[0], value));

		return *m;

	}

	
	/**
	 * @brief        Get reference to a complex single matrix
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	template <class T> Matrix<T>& 
	Get              (const string name) {

		reflist::iterator it = m_ref.find(name);
		return *boost::any_cast<Ptr<Matrix<T> > >(m_store[it->second[0]]);

	}

	
	/**
	 * @brief        Remove a complex double matrix
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	inline bool 
	Free             (const string name) {

		reflist::iterator nit = m_ref.find(name);
		
		if (nit == m_ref.end())
			return false;

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

		return true;
	}
	

 private:

	/** 
	 * @brief Private constructor (access Instance()). 
	 *        Singleton object for data storage. Access through Instance().
	 */
    Workspace () {}; 

	/**
	 * @brief Private copy constructor (access Instance()).
	 */
	Workspace (const Workspace&) {};

	reflist m_ref; /*! @brief Names and hash tags            */
	store   m_store;

	static Workspace *m_inst; /*!< @brief Single database instance       */
	
};

#endif /* _WORK_SPACE_H_ */
