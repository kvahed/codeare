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
	 * @brief        Add a complex matrix to complex matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	Matrix<cxfl>& 
	AddCXFL          (const string name, Ptr< Matrix<cxfl> > m);
	
	
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
	 * @brief        Get data from recon (Remote access)
	 *
	 * @param  name  Name
	 * @param  c     Raw data storage 
	 */
	void
	GetCXFL          (const string name, cxfl_data& c);
	
	
	/**
	 * @brief        Set data for recon (Remote access)
	 *
	 * @param name   Name
	 * @param c      CXFL data
	 */
	void 
	SetCXFL          (const string name, const cxfl_data& c);
	
	
	/*
	 * @brief        Add a complex matrix to complex matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	Matrix<cxdb>& 
	AddCXDB          (const string name, Ptr< Matrix<cxdb> > m);
	
	
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
	 * @brief        Get data from recon (Remote access)
	 *
	 * @param  name  Name
	 * @param  c     Raw data storage 
	 */
	void
	GetCXDB          (const string name, cxdb_data& c);
	
	
	/**
	 * @brief        Set data for recon (Remote access)
	 *
	 * @param name   Name
	 * @param c      CXDB data
	 */
	void 
	SetCXDB          (const string name, const cxdb_data& c);
	
	
	/**
	 * @brief        Add a real matrix to real matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Reference to 
	 */
	Matrix<float>&
	AddRLFL          (const string name, Ptr< Matrix<float> > m) ;
	
	
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
	 * @brief        Get data from recon (Remote access)
	 *
	 * @param  name  Name
	 * @param  r     RLDB data storage
	 */
	void 
	GetRLFL          (const string name, rlfl_data& r);
		

	/**
	 * @brief        Set data for recon (Local access)
	 * 
	 * @param  name  Name
	 * @param  r     RLDB data
	 */
	void 
	SetRLFL          (const string name, const rlfl_data& r);
		
	
	/**
	 * @brief        Add a real matrix to real matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Reference to 
	 */
	Matrix<double>&
	AddRLDB          (const string name, Ptr< Matrix<double> > m) ;
	
	
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
	 * @brief        Get data from recon (Remote access)
	 *
	 * @param  name  Name
	 * @param  r     RLDB data storage
	 */
	void 
	GetRLDB          (const string name, rldb_data& r);
		

	/**
	 * @brief        Set data for recon (Local access)
	 * 
	 * @param  name  Name
	 * @param  r     RLDB data
	 */
	void 
	SetRLDB          (const string name, const rldb_data& r);
		
	
	/**
	 * @brief        Add a short matrix to short matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	Matrix<short>&
	AddSHRT          (const string name, Ptr< Matrix<short> > m);

	
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
	 * @brief        Get data from recon
	 *
	 * @param  name  Name
	 * @param  p     Short data storage
	 */
	void 
	GetSHRT          (const string name, shrt_data& p);
	
	
	/**
	 * @brief        Set data for recon
	 *
	 * @param  name  Name
	 * @param  p     Short data
	 */
	void 
	SetSHRT          (const string name, const shrt_data& p);
	
	
	/**
	 * @brief        Add a long matrix to long matrix container
	 *
	 * @param  name  Name
	 * @param  m     The added matrix
	 * @return       Success
	 */
	Matrix<long>&
	AddLONG          (const string name, Ptr< Matrix<long> > m);

	
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
	 * @brief        Get data from recon
	 *
	 * @param  name  Name
	 * @param  p     Long data storage
	 */
	void 
	GetLONG          (const string name, long_data& p);
	
	
	/**
	 * @brief        Set data for recon
	 *
	 * @param  name  Name
	 * @param  p     Long data
	 */
	void 
	SetLONG          (const string name, const long_data& p);
	
	
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
