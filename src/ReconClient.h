/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#ifndef __RECON_CLIENT_H__
#define __RECON_CLIENT_H__

#include <complex>
#include "Matrix.h"
#include "Configurable.h"


#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

using namespace RRSModule;

/**
 * @brief Remote recon client
 */
namespace RRClient {
	
	/**
	 * @brief              CORBA ICE reconstruction client 
	 */
	class ReconClient :
	public Configurable {
		
		
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
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Process             (const char* name);
		
		/**
		 * @brief           Set raw my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		SetRaw              (Matrix< std::complex<float> >& M) {
			
			m_raw->dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				m_raw->dims[j] = M.Dim(j);
			
			m_raw->dreal.length(M.Size()); 
			m_raw->dimag.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++) {
				m_raw->dreal[i] = M[i].real();
				m_raw->dimag[i] = M[i].imag(); 
			}
			
			m_rrsi->raw(*(m_raw));
			
		};
		
		
		/**
		 * @brief           Return raw data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		GetRaw              (Matrix< std::complex<float> >& M) {
			
			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = m_raw->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetRawSize(); i++)
				M[i] = std::complex<float>(m_raw->dreal[i],m_raw->dimag[i]);
			
		};
		
		
		/**
		 * @brief           Set helper my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void
		SetHelper              (Matrix< double >& M) {
			
			for (int j = 0; j < INVALID_DIM; j++)
				m_helper->dims[j] = M.Dim(j);
			
			m_helper->vals.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++)
				m_helper->vals[i] = M[i];
			
			m_rrsi->helper(*(m_helper));
			
		};
		
		
		/**
		 * @brief           Return helper data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void
		GetHelper              (Matrix< double >& M) {
			
			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = m_helper->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetHelperSize(); i++)
				M[i] = m_helper->vals[i];
			
		};
		
		
		/**
		 * @brief           Set my Pixel data
		 * 
		 * @param  M        Given matrix
		 */
		void 
		SetPixel            (Matrix<short>& M) {
			
			for (int j = 0; j < INVALID_DIM; j++)
				m_pixel->dims[j] = M.Dim(j);
			
			m_pixel->vals.length(M.Size()); 
			
			for (int i = 0; i < M.Size(); i++)
				m_pixel->vals[i] = M[i];
			
			m_rrsi->pixel(*(m_pixel));
			
		};
		
		
		/**
		 * @brief           Put my Pixel data into 
		 * 
		 * @param  M        Given matrix
		 */
		void 
		GetPixel            (Matrix<short>& M) {
			
			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = m_pixel->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetPixelSize(); i++)
				M[i] = m_pixel->vals[i];
			
		};
		
		
		/**
		 * @brief           Raw repository size
		 *
		 * @return          Size
		 */
		long
		GetRawSize             ();
		
		
		/**
		 * @brief           Pixel repository size
		 *
		 * @return          Size
		 */
		long
		GetPixelSize           ();
		
		
		/**
		 * @brief           Helper repository size
		 *
		 * @return          Size
		 */
		long
		GetHelperSize          ();
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< Remote Recon interface               */
		
		raw_data*           m_raw;        /**< Raw data    (complex float sequence) */
		helper_data*        m_helper;     /**< Helper data (complex float sequence) */
		pixel_data*         m_pixel;      /**< Pixel data  (short sequence)         */
		
		CORBA::ORB_var      m_orb;        /**< Orb                                  */
		
	};
	
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
		std::cerr << "CORBA COMM_FAILURE" << std::endl; 
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
		std::cerr << "CORBA Exception" << std::endl; 
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
		std::cerr << "CORBA Fatal Exception" << std::endl; 
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
		std::cerr << "Exception" << std::endl; 
	};
};

#endif // __RECON_CLIENT_H__
