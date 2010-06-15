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
	 * @brief             Construct and initialise remote recon interface
	 *                    CORBA connection is initialised, data sequences setup.
	 *                    Set of exception handling for different stages of errors.
	 *                    Thus far only omniORB4 supported.
	 *          
	 * @param  name       Service name announced to name server
	 * @param  tracelevel Trace level [0-40]
	 *                    level 0:  critical errors only
     *                    level 1:  informational messages only
	 *					  level 2:  configuration information and warnings
	 *					  level 5:  notifications when server threads are created and communication endpoints are shutdown
	 *					  level 10: execution and exception traces
	 *					  level 25: trace each send or receive of a giop message
	 *					  level 30: dump up to 128 bytes of each giop message
	 *					  level 40: dump complete contents of each giop message
	 */
	ReconClient           (const char* name, const char* tracelevel);
	
	
	/**
	 * @brief             Destruct CORBA data and connection 
	 */
	~ReconClient          ();
	

	/**
	 * @brief             Request server-side data procession 
	 *
	 * @param  method     Recon method
	 *
	 * @return            Error code
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


	/**
	 * @brief           Set my pixel data with ...
	 *
	 * @param  M        Given matrix
	 */
	void 
	SetPixel            (Matrix<short>& M) {
		
		shorts vals;
		int    i;
		
		m_pixel->dims.length(16);
		
		for (i = 0; i < 16; i++)
			m_pixel->dims[i] = M.Dim(i);
		
		long size = GetPixelSize();
		
		m_pixel->vals.length(size); 
		
		for (i = 0; i < size; i++)
			m_pixel->vals[i] = M[i];
		
		m_rrsi->pixel(m_pixel[0]);
	};
	

	/**
	 * @brief           Return pixel data after recon to ...
	 *
	 * @param  M        Given matrix
	 */
	void 
	GetPixel            (Matrix<short>& M) {
	
		short* pixel = new short [GetPixelSize()];
		int             dim[INVALID_DIM], i;

		for (i = 0; i < INVALID_DIM; i++)
			dim[i] = m_pixel->dims[i];
		
		M.Reset(pixel, dim);
		for (i = 0; i < GetPixelSize(); i++)
			M[i] = (short)m_pixel->vals[i];
	
	};
	

	/**
	 * @brief               Raw repository size
	 *
	 * @return              Size
	 */
	long
	GetRawSize              ();
	

	/**
	 * @brief               Pixel repository size
	 *
	 * @return              Size
	 */
	long
	GetPixelSize            ();
	

	/**
	 * @brief               Helper repository size
	 *
	 * @return              Size
	 */
	long
	GetHelperSize           ();
	
	

private:
	

	RRSInterface_var             m_rrsi;   /**< My remote recon interface                       */
	
	raw_data*                    m_raw;    /**< Raw data repository   (complex float sequence)  */
	raw_data*                    m_helper; /**< Helper data repository (complex float sequence) */
	pixel_data*                  m_pixel;  /**< Pixel data repository (short sequence)          */

	CORBA::ORB_var               m_orb;    /**< Object request broker                           */
	
};


/**
 * @brief CORBA server connection exception
 */
class DS_ServerConnectionException {
public:
	/**
	 * @brief Handle CORBA server connection exceptions
	 */
	DS_ServerConnectionException () { 
		cerr << "CORBA COMM_FAILURE" << endl; 
	};
};


/**
 * @brief CORBA system exception
 */
class DS_SystemException           {
public:
	/**
	 * @brief Handle system exceptions
	 */
	DS_SystemException           () { 
		cerr << "CORBA Exception" << endl; 
	};
};


/**
 * @brief CORBA fatal exception
 */
class DS_FatalException            {
public:
	/**
	 * @brief Handle fatals exceptions
	 */
	DS_FatalException            () { 
		cerr << "CORBA Fatal Exception" << endl; 
	};
};


/**
 * @brief CORBA general exception :(
 */
class DS_Exception                 {
public:
	/**
	 * @brief Not much help here
	 */
	DS_Exception                 () { 
		cerr << "Exception" << endl; 
	};
};

#endif // __RECON_CLIENT_H__
