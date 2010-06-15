#ifndef __RECON_CLIENT_H__
#define __RECON_CLIENT_H__

#include <complex>

#include "Matrix.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

using namespace std;
using namespace RRSModule;
                                                                                

/**
 * @brief              CORBA ICE reconstruction client 
 */
class ReconClient                  {
	
	
public:
	
	/**
	 * @brief           Construct and initialise remote recon interface
	 */
	ReconClient         (const char* name, const char* tracelevel);
	
	
	/**
	 * @brief           Default destructor
	 */
	~ReconClient        ();
	
	/**
	 * @brief           Request data procession on recon service
	 *
	 * @param  method   Recon method
	 * @return          Error code
	 */ 
	error_code              
	Requestprocess_data (int method);
	
	/**
	 * @brief           Set raw my data with ...
	 *
	 * @param  M        Given matrix
	 */
	void 
	SetRaw              (Matrix< complex<float> >& M) {
		
		floats dabs;
		floats darg;
		int    i;
		
		for (i = 0; i < 16; i++)
			m_raw->dims[i] = M.Dim()[i];
		
		long size = GetSize();
		
		cout << size << endl;
		
		m_raw->dabs.length(size); 
		m_raw->darg.length(size);
		
		for (i = 0; i < size; i++) {
			m_raw->dabs[i] = abs(M[i]);
			m_raw->darg[i] = arg(M[i]); 
		}
		
		m_rrsi->raw(m_raw[0]);
		
	};
	

	/**
	 * @brief           Return raw data after recon to ...
	 *
	 * @param  M        Given matrix
	 */
	void 
	GetRaw              (Matrix< complex<float> >& M) {
		
		complex<float>* raw = new complex<float> [GetSize()];
		int             dim[INVALID_DIM], i;

		for (i = 0; i < INVALID_DIM; i++)
			dim[i] = m_raw->dims[i];
		
		M.Reset(raw, dim);
		for (i = 0; i < GetSize(); i++)
			M[i] = complex<float>(m_raw->dabs[i],m_raw->darg[i]);
		
	};
	

	void 
	SetPixel            (Matrix<short>* M) {};
	
	void 
	GetPixel            (Matrix<short>* M) {};
	
	/**
	 * @brief Repository size
	 *
	 * @return          Size
	 */
	long
	GetSize             ();
	
	
private:
	
	RRSInterface_var             m_rrsi;   /**< Remote Recon interface               */
	
	raw_data*                    m_raw;    /**< Raw data    (complex float sequence) */
	raw_data*                    m_helper; /**< Helper data (complex float sequence) */
	pixel_data*                  m_pixel;  /**< Pixel data  (short sequence)         */

	CORBA::ORB_var               m_orb;
	
};


class DS_ServerConnectionException {
public:
	DS_ServerConnectionException () { 
		cerr << "CORBA COMM_FAILURE" << endl; 
	};
};


class DS_SystemException           {
public:
	DS_SystemException           () { 
		cerr << "CORBA Exception" << endl; 
	};
};


class DS_FatalException            {
public:
	DS_FatalException            () { 
		cerr << "CORBA Fatal Exception" << endl; 
	};
};


class DS_Exception                 {
public:
	DS_Exception                 () { 
		cerr << "Exception" << endl; 
	};
};

#endif // __RECON_CLIENT_H__
