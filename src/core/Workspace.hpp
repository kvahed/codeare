#ifndef __WORK_SPACE_HPP_
#define __WORK_SPACE_HPP_

#include "DllExport.h"
#include "Matrix.hpp"
#include "Configurable.hpp"
#include "Params.hpp"

#include <boost/any.hpp>
#ifdef HAVE_CXX11_SHARED_PTR
#  include <memory>
#  define  shrd_ptr std::shared_ptr
#  define  mk_shared std::make_shared
#else
#  include <boost/shared_ptr.hpp>
#  include <boost/make_shared.hpp>
#  define  shrd_ptr boost::shared_ptr
#  define  mk_shared boost::make_shared
#endif

#include <map>

#ifdef __APPLE__
  #include "AppleDigest.hpp"
#else
  #include "Digest.hpp"
#endif

#include <boost/property_tree/ptree.hpp>

typedef map<string, std::vector<std::string> > reflist;
typedef pair<string, std::vector<std::string> > refent;
typedef map<string, boost::any> store;
typedef pair<string, boost::any> entry;


template<class T> struct PrintTraits;


/**
 * @brief Global workspace. Singleton.
 */
class DLLEXPORT Workspace {



 public:


	/**
	 * @brief        Clean up
	 */
	virtual ~Workspace        ();


	/**
	 * @brief        Get reference to database instance
	 */
	static Workspace& Instance  () {
        if (m_inst == 0)
            m_inst = new Workspace ();        
        return *m_inst;
    }


	/**
	 * @brief        Initialise database
	 */
	codeare::error_code 
	Initialise       ();


	/**
	 * @brief        Free RAM
	 *
	 * @return       Success
	 */ 
	codeare::error_code 
	Finalise         ();
	
	
    template <class T> inline Matrix<T>&
	Get              (const std::string& name) throw () {

		reflist::iterator it = m_ref.find(name);

        if (it == m_ref.end())
            printf ("*** WARNING: Matrix %s could not be found in workspace!\n", name.c_str());

        const boost::any& ba = m_store[it->second[0]];

        try {
			boost::any_cast<shrd_ptr<Matrix<T> > >(m_store[it->second[0]]);
		} catch (const boost::bad_any_cast& e) {
			printf ("*** WARNING: Failed to retrieve %s - %s.\n             Requested %s - have %s.\n",
					name.c_str(), e.what(),
                    demangle(typeid(shrd_ptr<Matrix<T> >).name()).c_str(),
                    demangle(ba.type().name()).c_str());
		}

		return *boost::any_cast<shrd_ptr<Matrix<T> > >(m_store[it->second[0]]);

	}

	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> inline codeare::error_code
	GetMatrix          (const std::string& name, Matrix<T>& m) {
		codeare::error_code ec = Exists<T>(name);
		if (ec == codeare::OK)
        	m = Get<T>(name);
		return ec;

	}

	
	
	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> inline void
	SetMatrix          (const std::string& name, Matrix<T>& m) {

	    std::vector<std::string> tag(2);
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();

		shrd_ptr<Matrix<T> > pm = mk_shared<Matrix<T> >();
		boost::any val     = pm;

		if (m_ref.find (name) == m_ref.end()) {
			m_ref.insert (refent(name,tag));
			m_store.insert (entry(tag[0], val));
		} else 

		m.SetClassName(name.c_str());

		*pm = m;

	}
    template<class T> inline void
    Set (const std::string& name, Matrix<T>& m) {
        SetMatrix (name, m);
    }
    template<class T> inline void
    Add (const std::string& name, Matrix<T>& m) {
        SetMatrix (name, m);
    }
	
	
	/**
	 * @brief        Add a matrix to workspace
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
     *
	 * @return       Success
	 */
	template<class T> inline Matrix<T>&
	AddMatrix        (const std::string& name, shrd_ptr< Matrix<T> > m) {

	  std::vector<std::string> tag(2);
		boost::any value = m;
        
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();
		reflist::iterator ri = m_ref.find (name);
		if (ri != m_ref.end())
			Free (name);
		m_ref.insert (refent(name, tag));
		m_store.insert (entry (tag[0], value));

        m->SetClassName(name.c_str());
        
		return *m;

	}

	/**
	 * @brief        Add a matrix to workspace
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
     *
	 * @return       Success
	 */
	template<class T> inline Matrix<T>&
	AddMatrix        (const std::string& name) {

		shrd_ptr<Matrix<T> > m = mk_shared<Matrix<T> >();
		std::vector<std::string> tag(2);
		boost::any value = m;
        
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();
		reflist::iterator ri = m_ref.find (name);
		if (ri != m_ref.end())
			Free (name);
		m_ref.insert (refent(name, tag));
		m_store.insert (entry (tag[0], value));

        m->SetClassName(name.c_str());
        
		return *m;

	}


	/**
	 * @brief        Add a matrix to workspace
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
     *
	 * @return       Success
	 */
	template<class T> inline Matrix<T>&
	AddMatrix        (const std::string& name, Matrix<T>& m) {
		return AddMatrix(name, &m);
	}
	
	/**
	 * @brief        Remove a complex double matrix
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	template<class T> inline codeare::error_code
	Exists (const std::string& name) const {

        reflist::const_iterator nit = m_ref.find(name);
        if (nit == m_ref.end())
        	return codeare::NO_MATRIX_IN_WORKSPACE_BY_NAME;
        else if (nit->second[1].compare(typeid(T).name()))
        	return codeare::WRONG_MATRIX_TYPE;
        return codeare::OK;

    }


 	/**
	 * @brief        Remove a complex double matrix
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	inline bool
	Free (const std::string& name) {
        
        reflist::iterator nit = m_ref.find(name);
        
        if (nit == m_ref.end())
            return false;
        
        store::iterator dit = m_store.find (nit->second[0]);
        
        if      (nit->second[1].compare(typeid(cxfl).name())   == 0)
            boost::any_cast<shrd_ptr<Matrix<cxfl  > > >(dit->second).reset();
        else if (nit->second[1].compare(typeid(cxdb).name())   == 0)
            boost::any_cast<shrd_ptr<Matrix<cxdb  > > >(dit->second).reset();
        else if (nit->second[1].compare(typeid(float).name())  == 0)
            boost::any_cast<shrd_ptr<Matrix<float > > >(dit->second).reset();
        else if (nit->second[1].compare(typeid(double).name()) == 0)
            boost::any_cast<shrd_ptr<Matrix<double> > >(dit->second).reset();
        else if (nit->second[1].compare(typeid(short).name())   == 0)
            boost::any_cast<shrd_ptr<Matrix<short > > >(dit->second).reset();
        else if (nit->second[1].compare(typeid(long).name())   == 0)
            boost::any_cast<shrd_ptr<Matrix<long  > > >(dit->second).reset();
        
        m_store.erase(dit);
        m_ref.erase(nit);
        
        return true;
        
    }
    

    /**
     * @brief        Get string representation of mapping
     *
     * @return       String representation of workspace content
     */
    void
    Print           (std::ostream& os) const;
    
    /**
     * @brief        Get casted parameter
     *
     * @param  key   Key
     * @return       Parameter
     */
    template<class T> inline T
    PGet (const std::string& key) const {
    	return p.Get<T>(key);
    }


    /**
	 * @brief        Get casted parameter
	 *
	 * @param  key   Key
	 * @return       Parameter
	 */
    template<class T> inline T
    PGet (const char* key) const {
        return p.Get<T>(std::string(key));
    }


    /**
	 * @brief        Set parameter
	 *
	 * @param  key   Key
	 * @param  val   Value
	 */
    inline void
    PSet (const std::string& key, const boost::any& val) {
    	p.Set(key,val);
    }


    /**
	 * @brief        Set parameter
	 *
	 * @param  key   Key
	 * @param  val   Value
	 */
    inline void
    PSet (const char* key, const boost::any& val) {
        p.Set(std::string(key),val);
    }


	Params p;        /** < @brief Global parameter list */

 private:

	/** 
	 * @brief        Private constructor (access Instance()). 
	 *               Singleton object for data storage. Access through Instance().
	 */
    Workspace        ();
    
	/**
	 * @brief        Private copy constructor (access Instance()).
	 */
	Workspace        (const Workspace&) {};

#ifdef _MSC_VER
#pragma warning (disable : 4251)
#endif
  reflist m_ref;   /**< @brief Names and hash tags               */
	store   m_store; /**< @brief Data pointers                     */
#ifdef _MSC_VER
#pragma warning (default : 4251)
#endif

	static Workspace* m_inst; /**< @brief Single database instance */

};

static Workspace& wspace = Workspace::Instance();


/**
 * @brief            Dump to ostream
 *
 * @param  os        Output stream
 * @param  w         Workspace
 * @return           The output stream
 */
inline static std::ostream&
operator<< (std::ostream& os, const Workspace& w) {
	w.Print(os);
    return os;
}



#endif /* _WORK_SPACE_H_ */
