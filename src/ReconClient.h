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
	 * @brief           Default destructor
	 */
	void 
	Cleanup             ();
	
	/**
	 * @brief           Request data procession on recon service
	 *
	 * @param  m        Recon method
	 * @return          Error code
	 */ 
	error_code              
	Requestprocess_data (method m);
	
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

		m_have_raw = true;
		
		m_raw->dims.length(INVALID_DIM);
		
		for (i = 0; i < INVALID_DIM; i++)
			m_raw->dims[i] = M.Dim(i);

		long size = GetRawSize();

		m_raw->dabs.length(size); 
		m_raw->darg.length(size);
		
		for (i = 0; i < size; i++) {
			m_raw->dabs[i] = M[i].real();
			m_raw->darg[i] = M[i].imag(); 
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
		
		int             dim[INVALID_DIM], i;

		for (i = 0; i < INVALID_DIM; i++)
			dim[i] = m_raw->dims[i];
		
		M.Reset(dim);
		for (i = 0; i < GetRawSize(); i++)
			M[i] = complex<float>(m_raw->dabs[i],m_raw->darg[i]);
		
	};
	

	/**
	 * @brief           Set helper my data with ...
	 *
	 * @param  M        Given matrix
	 */
	void
	SetHelper              (Matrix< complex<float> >& M) {

		floats dabs;
		floats darg;
		int    i;

		m_have_helper = true;

		m_helper->dims.length(INVALID_DIM);

		for (i = 0; i < INVALID_DIM; i++)
			m_helper->dims[i] = M.Dim(i);

		long size = GetHelperSize();

		m_helper->dabs.length(size);
		m_helper->darg.length(size);

		for (i = 0; i < size; i++) {
			m_helper->dabs[i] = M[i].real();
			m_helper->darg[i] = M[i].imag();
		}

		m_rrsi->helper(m_helper[0]);

	};


	/**
	 * @brief           Return helper data after recon to ...
	 *
	 * @param  M        Given matrix
	 */
	void
	GetHelper              (Matrix< complex<float> >& M) {

		int             dim[INVALID_DIM], i;

		for (i = 0; i < INVALID_DIM; i++)
			dim[i] = m_helper->dims[i];

		M.Reset(dim);
		for (i = 0; i < GetHelperSize(); i++)
			M[i] = complex<float>(m_helper->dabs[i],m_helper->darg[i]);

	};


	/**
	 * @brief           Set my Pixel data
	 * 
	 * @param  M        Given matrix
	 */
	void 
	SetPixel            (Matrix<short>& M) {
		
		shorts vals;
		int    i;
		
		m_have_pixel = true;
		
		m_pixel->dims.length(INVALID_DIM);
		
		for (i = 0; i < INVALID_DIM; i++)
			m_pixel->dims[i] = M.Dim(i);
		
		long   size = GetPixelSize();
		
		m_pixel->vals.length(size); 
		
		for (i = 0; i < size; i++)
			m_pixel->vals[i] = M[i];
		
		m_rrsi->pixel(m_pixel[0]);
	};
	

	/**
	 * @brief           Put my Pixel data into 
	 * 
	 * @param  M        Given matrix
	 */
	void 
	GetPixel            (Matrix<short>& M) {
	
		int             dim[INVALID_DIM], i;

		for (i = 0; i < INVALID_DIM; i++)
			dim[i] = m_pixel->dims[i];
		
		M.Reset(dim);
		for (i = 0; i < GetPixelSize(); i++)
			M[i] = (short)m_pixel->vals[i];
	
	};
	
	/**
	 * @brief Repository size
	 *
	 * @return          Size
	 */
	long
	GetRawSize             ();
	
	/**
	 * @brief Repository size
	 *
	 * @return          Size
	 */
	long
	GetPixelSize             ();
	
	/**
	 * @brief Repository size
	 *
	 * @return          Size
	 */
	long
	GetHelperSize             ();
	
	
private:
	
	RRSInterface_var             m_rrsi;   /**< Remote Recon interface               */
	
	raw_data*                    m_raw;    /**< Raw data    (complex float sequence) */
	raw_data*                    m_helper; /**< Helper data (complex float sequence) */
	pixel_data*                  m_pixel;  /**< Pixel data  (short sequence)         */

	bool                         m_have_raw;
	bool                         m_have_pixel;
	bool                         m_have_helper;

	strings*                     m_labels;

	CORBA::ORB_var               m_orb;
	
};


/**
 * @brief Corba communication exception
 */
class DS_ServerConnectionException {
public:
	/**
	 * @brief Handle communication exceptions
	 */
	DS_ServerConnectionException () { 
		cerr << "CORBA COMM_FAILURE" << endl; 
	};
};


/**
 * @brief Corba system exception
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
 * @brief Corba fatal exception
 */
class DS_FatalException            {
public:
	/**
	 * @brief Handle fatal exception
	 */
	DS_FatalException            () { 
		cerr << "CORBA Fatal Exception" << endl; 
	};
};


/**
 * @brief Corba unspecified exception
 */
class DS_Exception                 {
public:
	/**
	 * @brief Handle unspecified exception
	 */
	DS_Exception                 () { 
		cerr << "Exception" << endl; 
	};
};

#endif // __RECON_CLIENT_H__
