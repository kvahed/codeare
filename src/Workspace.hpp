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

template<class T>
struct CorbaTraits {};

template<>
struct CorbaTraits<cxfl> {
	typedef RRSModule::cxfl_data Type;
};
template<>
struct CorbaTraits<cxdb> {
	typedef RRSModule::cxdb_data Type;
};
template<>
struct CorbaTraits<float> {
	typedef RRSModule::rlfl_data Type;
};
template<>
struct CorbaTraits<double> {
	typedef RRSModule::rldb_data Type;
};
template<>
struct CorbaTraits<long> {
	typedef RRSModule::long_data Type;
};
template<>
struct CorbaTraits<short> {
	typedef RRSModule::shrt_data Type;
};

/**
 * @brief Central database for all shared matrices<br/>
 *        Matrix smart pointers are stored in individual (type) maps along with string identifier and retrieved by their names.
 */
class Workspace : public Configurable {


	typedef map<string, string[2]> reflist;
	typedef pair<string, string[2]> ref;
	typedef map<string, Ptr< Matrix<cxfl> > > cxfl_db;
	typedef map<string, Ptr< Matrix<cxdb> > > cxdb_db;
	typedef map<string, Ptr< Matrix<float> > > rlfl_db;
	typedef map<string, Ptr< Matrix<double> > > rldb_db;
	typedef map<string, Ptr< Matrix<short> > > shrt_db;
	typedef map<string, Ptr< Matrix<long> > > long_db;
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
	 * @brief        Get data from recon (Remote connector)
	 *
	 * @param  name  Name
	 * @param  t     Raw data storage type
	 */
	template <class T> void
	GetMatrix        (const string& name, typename CorbaTraits<T>::Type& t) {

		if (m_ref.find (name) == m_ref.end())
			return;

		map<string,string[2]>::iterator it = m_ref.find(name);

		Ptr< Matrix<cxfl> > tmp = boost::any_cast<Ptr<Matrix<cxfl> > >(m_store[it->second[0]]);

		for (int j = 0; j < INVALID_DIM; j++) {
			t.dims[j] = tmp->Dim(j);
			t.res[j]  = tmp->Res(j);
		}

		t.vals.length(2 * tmp->Size());

		memcpy (&t.vals[0], &tmp->At(0), tmp->Size() * sizeof(T));

	}
	
	
	/**
	 * @brief        Set data for recon (Remote connector)
	 *
	 * @param name   Name
	 * @param t      Raw data storage type
	 */
	template <class T> void 
	SetMatrix        (const string name, const typename CorbaTraits<T>::Type& t) {

		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = t.dims[i];
			mress[i] = t.res[i];
		}

		Ptr<Matrix<T> > pm = NEW (Matrix<T>(mdims, mress));
		boost::any val = pm;

		if (m_ref.find (name) != m_ref.end())
			Free (name);

		std::string tag[2];
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();

		m_ref.insert (pair<string, string[2]> (name, tag));
		m_store.insert (entry (tag[0], val));

		pm->SetClassName(name.c_str());

		memcpy (&pm->At(0), &t.vals[0], pm->Size() * sizeof(T));

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

		map<string,string[2]>::iterator it = m_ref.find(name);
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
/*
		reflist::iterator nit = m_ref.find(name);
		
		if (nit == m_ref.end())
			return false;


		if        (nit->second[1].compare("cxfl") == 0) {
			store::iterator dit = m_store.find (nit->second[0]);
			delete boost::any_cast<Ptr<Matrix<cxfl> > >(dit->second);
			m_store.erase(dit);
		} else if (nit->second[1].compare("cxdb") == 0) {
			cxdb_db::iterator dit = m_cxdb.find (nit->second[0]);
			delete dit->second;
			m_cxdb.erase(dit);
		} else if (nit->second[1].compare("float") == 0) {
			rlfl_db::iterator dit = m_rlfl.find (nit->second[0]);
			delete dit->second;
			m_rlfl.erase(dit);
		} else if (nit->second[1].compare("double") == 0) {
			rldb_db::iterator dit = m_rldb.find (nit->second[0]);
			delete dit->second;
			m_rldb.erase(dit);
		} else if (nit->second[1].compare("shrt") == 0) {
			shrt_db::iterator dit = m_shrt.find (nit->second[0]);
			delete dit->second;
			m_shrt.erase(dit);
		} else if (nit->second[1].compare("long") == 0) {
			long_db::iterator dit = m_long.find (nit->second[0]);
			delete dit->second;
			m_long.erase(dit);
		} 

		m_ref.erase(nit);
		*/
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
