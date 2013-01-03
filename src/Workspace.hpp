#ifndef __WORK_SPACE_HPP_
#define __WORK_SPACE_HPP_

#include "Matrix.hpp"
#include "Configurable.hpp"

#include <map>

#ifdef __WIN32__ 
  #include "RRSModule.h"
#else
  #include "RRSModule.hh"
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
	typedef map<string, Ptr< Matrix<cxfl> > > cxfl_db;
	typedef map<string, Ptr< Matrix<cxdb> > > cxdb_db;
	typedef map<string, Ptr< Matrix<float> > > rlfl_db;
	typedef map<string, Ptr< Matrix<double> > > rldb_db;
	typedef map<string, Ptr< Matrix<short> > > shrt_db;
	typedef map<string, Ptr< Matrix<long> > > long_db;

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
	GetMatrix          (const string name, Matrix<T>& m);
	
	
	/**
	 * @brief        Get data from recon (Local connector)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> void
	SetMatrix          (const string name, Matrix<T>& m);
	
	
	/**
	 * @brief        Get data from recon (Remote connector)
	 *
	 * @param  name  Name
	 * @param  t     Raw data storage type
	 */
	template <class T> void
	GetMatrix        (const string name, T& t);
	
	
	/**
	 * @brief        Set data for recon (Remote connector)
	 *
	 * @param name   Name
	 * @param t      Raw data storage type
	 */
	template <class T> void 
	SetMatrix        (const string name, const T& t);
	
	
	/**
	 * @brief        Add a matrix to according container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	template <class T> Matrix<T>& 
	AddMatrix        (const string name, Ptr< Matrix<T> > m);
	
	
	/***
	 * @brief        Add a matrix to according container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	//template <class T> Matrix<T>& 
	//AddMatrix        (const string& name, const data_type dt, Ptr< Matrix<T> > m);
	
	
	/**
	 * @brief        Get reference to a complex single matrix
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	template <class T> Matrix<T>& 
	Get              (const string name);
	
	
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

		if        (nit->second[1].compare("cxfl") == 0) {
			cxfl_db::iterator dit = m_cxfl.find (nit->second[0]);
			delete dit->second;
			m_cxfl.erase(dit);
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
		
		return true;
	}
	
	
	/**
	 * @brief        Get reference to complex single store
	 *
	 * @return       Reference to complex single store
	 */
	map < string, Ptr< Matrix<cxfl> > >& 
	CXFLMap          ();
		

	/**
	 * @brief        Get reference to complex single store
	 *
	 * @return       Reference to complex single store
	 */
	map < string, Ptr< Matrix<cxdb> > >& 
	CXDBMap          ();
		

	/**
	 * @brief        Get reference to complex double store
	 *
	 * @return       Reference to complex double store
	 */
	map < string, Ptr< Matrix<float> > >& 
	RLFLMap          ();
		

	/**
	 * @brief        Get reference to real single store
	 *
	 * @return       Reference to real single store
	 */
	map < string, Ptr< Matrix<double> > >& 
	RLDBMap          ();
		

	/**
	 * @brief        Get reference to short int store
	 *
	 * @return       Reference to short int store
	 */
	map < string, Ptr< Matrix<short> > >& 
	SHRTMap          ();


	/**
	 * @brief        Get reference to long int store
	 *
	 * @return       Reference to long int store
	 */
	map < string, Ptr< Matrix<long> > >& 
	LONGMap          ();
		

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

	reflist             m_ref; /*! @brief Names and hash tags            */
	
	cxfl_db m_cxfl; /*!< @brief Complex float data repository  */
	cxdb_db m_cxdb; /*!< @brief Complex double data repository */
	rlfl_db m_rlfl; /*!< @brief Real float data repository     */
	rldb_db m_rldb; /*!< @brief Real double data repository    */
	shrt_db m_shrt; /*!< @brief Integer data respository       */
	long_db m_long; /*!< @brief Integer data respository       */
	
	static Workspace *m_inst; /*!< @brief Single database instance       */
	
};

#endif /* _WORK_SPACE_H_ */
