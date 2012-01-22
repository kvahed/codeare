#ifndef _DATA_BASE_H_
#define _DATA_BASE_H_

#include "Matrix.hpp"

#ifdef __WIN32__ 
  #include "RRSModule.h"
#else
  #include "RRSModule.hh"
#endif

using namespace std;
using namespace RRSModule;

/**
 * @brief Central database for all shared matrices
 */
class DataBase {



 public:


	/**
	 * @brief        Clean up
	 */
	virtual ~DataBase        ();


	/**
	 * @brief        Get pointer to database instance
	 */
	static DataBase* 
	Instance         () ;
	

	/**
	 * @brief        Initialise database
	 */
	const error_code 
	Initialise       ();


	/**
	 * @brief        Free RAM
	 *
	 * @return       Success
	 */ 
	const error_code 
	Finalise         ();
	
	
	/**
	 * @brief        Get data from recon (Local access)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> void
	GetMatrix          (const string name, Matrix<T>& m);
	
	
	/**
	 * @brief        Get data from recon (Local access)
	 *
	 * @param  name  Name
	 * @param  m     CXFL data storage 
	 */
	template <class T> void
	SetMatrix          (const string name, Matrix<T>& m);
	
	
	/**
	 * @brief        Get data from recon (Remote access)
	 *
	 * @param  name  Name
	 * @param  c     Raw data storage 
	 */
	template <class T> void
	GetMatrix        (const string name, T& t);
	
	
	/**
	 * @brief        Set data for recon (Remote access)
	 *
	 * @param name   Name
	 * @param c      CXFL data
	 */
	template <class T> void 
	SetMatrix        (const string name, const T& t);
	
	
	/**
	 * @brief        Add a complex matrix to complex matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	template <class T> Matrix<T>& 
	AddMatrix        (const string name, Ptr< Matrix<T> > m);
	
	
	/**
	 * @brief        Remove a complex matrix from complex container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeCXFL         (const string name);
	
	
	/**
	 * @brief        Get reference to complex data by name
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	Matrix<cxfl>& 
	GetCXFL          (const string name);
	
	
	/**
	 * @brief        Remove a complex matrix from complex container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeCXDB         (const string name);
	
	
	/**
	 * @brief        Get reference to complex data by name
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	Matrix<cxdb>& 
	GetCXDB          (const string name);
	
	
	/**
	 * @brief        Remove a real matrix from real container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeRLFL         (const string name);
		
	
	/**
	 * @brief        Get reference to complex data by name
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	Matrix<float>& 
	GetRLFL          (const string name);
		
		
	/**
	 * @brief        Remove a real matrix from real container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeRLDB         (const string name);
		
	
	/**
	 * @brief       Get reference to complex data by name
	 * 
	 * @param  name Name
	 * @return      Reference to data if existent
	 */
	Matrix<double>& 
	GetRLDB         (const string name);
		
		
	/**
	 * @brief        Remove a short matrix from short container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeSHRT         (const string name);
	

	/**
	 * @brief        Get reference to complex data by name
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	Matrix<short>&   
	GetSHRT          (const string name);
	

	/**
	 * @brief        Remove a long matrix from long container
	 *
	 * @param  name  Name
	 * @return       Success
	 */
	bool 
	FreeLONG         (const string name);
	

	/**
	 * @brief        Get reference to complex data by name
	 * 
	 * @param  name  Name
	 * @return       Reference to data if existent
	 */
	Matrix<long>&   
	GetLONG          (const string name);
	

	/**
	 * @brief        Reference to complex single store
	 *
	 * @return       Reference to complex single store
	 */
	map < string, Ptr< Matrix<cxfl> > >& 
	CXFLMap          ();
		

	/**
	 * @brief        Reference to complex single store
	 *
	 * @return       Reference to complex single store
	 */
	map < string, Ptr< Matrix<cxdb> > >& 
	CXDBMap          ();
		

	/**
	 * @brief        Reference to complex double store
	 *
	 * @return       Reference to complex double store
	 */
	map < string, Ptr< Matrix<float> > >& 
	RLFLMap          ();
		

	/**
	 * @brief        Reference to real single store
	 *
	 * @return       Reference to real single store
	 */
	map < string, Ptr< Matrix<double> > >& 
	RLDBMap          ();
		

	/**
	 * @brief        Reference to short int store
	 *
	 * @return       Reference to short int store
	 */
	map < string, Ptr< Matrix<short> > >& 
	SHRTMap          ();


	/**
	 * @brief        Reference to long int store
	 *
	 * @return       Reference to long int store
	 */
	map < string, Ptr< Matrix<long> > >& 
	LONGMap          ();
		
 private:

	/** 
	 * @brief Private constructor. Use Instance
	 */
    DataBase() {}; 

	map < string, Ptr< Matrix<cxfl>   > > m_cxfl; /*!< @brief Complex float data repository  */
	map < string, Ptr< Matrix<cxdb>   > > m_cxdb; /*!< @brief Complex double data repository */
	map < string, Ptr< Matrix<float>  > > m_rlfl; /*!< @brief Real float data repository     */
	map < string, Ptr< Matrix<double> > > m_rldb; /*!< @brief Real double data repository    */
	map < string, Ptr< Matrix<short>  > > m_shrt; /*!< @brief Integer data respository       */
	map < string, Ptr< Matrix<long>   > > m_long; /*!< @brief Integer data respository       */
	
	static DataBase *m_inst; /*!< @brief Single database instance       */
	
};


#endif /* _DATA_BASE_H_ */
