#ifndef __WORK_SPACE_HPP_
#define __WORK_SPACE_HPP_

#include "Matrix.hpp"
#include "Configurable.hpp"
#include "Params.hpp"

#include <boost/any.hpp>
#include <map>

#ifdef __APPLE__
  #include "AppleDigest.hpp"
#else
  #include "Digest.hpp"
#endif

using namespace std;



typedef map<string, string[2]> reflist;
typedef pair<string, string[2]> ref;
typedef map<string, boost::any> store;
typedef pair<string, boost::any> entry;


/**
 * @brief Global workspace. Singleton.
 */
class Workspace {



 public:


	/**
	 * @brief        Clean up
	 */
	virtual ~Workspace        ();


	/**
	 * @brief        Get reference to database instance
	 */
	static Workspace&
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
	
	
    template <class T> inline Matrix<T>&
	Get              (const string& name) {

        reflist::iterator it = m_ref.find(name);

        if (it == m_ref.end())
            printf ("**WARNING**: Matrix %s could not be found in workspace!\n");

        const boost::any& ba = m_store[it->second[0]];
        
        try {
			boost::any_cast<Ptr<Matrix<T> > >(m_store[it->second[0]]);
		} catch (const boost::bad_any_cast& e) {
			printf ("**WARNING**: Failed to retrieve %s - %s.\n             Requested %s - have %s.\n",
					name.c_str(), e.what(),
                    demangle(typeid(Ptr<Matrix<T> >).name()).c_str(),
                    demangle(ba.type().name()).c_str());
		}

		return *boost::any_cast<Ptr<Matrix<T> > >(m_store[it->second[0]]);

	}

	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> inline void
	GetMatrix          (const string& name, Matrix<T>& m) {
        m = Get<T>(name);
	}

	
	
	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> inline void
	SetMatrix          (const string& name, Matrix<T>& m) {

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
    template<class T> inline void
    Set (const string& name, Matrix<T>& m) {
        SetMatrix (name, m);
    }
    template<class T> inline void
    Add (const string& name, Matrix<T>& m) {
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
	AddMatrix        (const string& name, Ptr< Matrix<T> > m = NEW (Matrix<T>)) {

		std::string tag[2];
		boost::any value = m;
        
		tag[0] = sha256(name);
		tag[1] = typeid(T).name();
		assert (m_ref.find (name) == m_ref.end());
		m_ref.insert (ref(name, tag));
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
	AddMatrix        (const string& name, Matrix<T>& m) {

		return AddMatrix(name, &m);

	}
	
	/**
	 * @brief        Remove a complex double matrix
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	inline bool
	Exists (const std::string& name) const {

        reflist::const_iterator nit = m_ref.find(name);
        return (nit == m_ref.end()) ? false : true;

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
    

    /**
     * @brief        Get string representation of mapping
     *
     * @return       String representation of workspace content
     */
    const char*
    c_str            () const {
    
        std::string sb;
        
        for (reflist::const_iterator i = m_ref.begin(); i != m_ref.end(); i++)
            sb += i->first + "\t" + i->second[0] + "\t" + i->second[1] + "\n";
        
        return sb.c_str();
        
    }


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

	reflist m_ref;   /**< @brief Names and hash tags               */
	store   m_store; /**< @brief Data pointers                     */

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
	os << w.c_str();
    return os;
}

#endif /* _WORK_SPACE_H_ */
