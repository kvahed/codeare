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

#include "Matrix.hpp"
#include "Configurable.hpp"


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
		 * @brief           Request data procession on recon service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Init                (const char* name);
		
 		/**
		 * @brief           Request data procession on recon service
		 *
		 * @param  name     Recon method
		 * @return          Error code
		 */ 
		error_code              
		Finalise            (const char* name);
		
		/**
		 * @brief           Set raw my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		SetRaw              (Matrix< std::complex<float> >& M) {
			
			raw_data r; 

			r.dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				r.dims[j] = M.Dim(j);
			
			r.dreal.length(M.Size()); 
			r.dimag.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++) {
				r.dreal[i] = M[i].real();
				r.dimag[i] = M[i].imag(); 
			}
			
			m_rrsi->raw(r);
	
		};
		
		
		/**
		 * @brief           Return raw data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		GetRaw              (Matrix< std::complex<float> >& M) {
			
			raw_data* rp = m_rrsi->raw();

			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = rp->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetSize(rp->dims); i++)
				M[i] = std::complex<float>(rp->dreal[i],rp->dimag[i]);
			
		};
		
		
		/**
		 * @brief           Set raw my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		SetRHelper              (Matrix< std::complex<float> >& M) {
			
			raw_data r;

			r.dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				r.dims[j] = M.Dim(j);
			
			r.dreal.length(M.Size()); 
			r.dimag.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++) {
				r.dreal[i] = M[i].real();
				r.dimag[i] = M[i].imag(); 
			}
			
			m_rrsi->rhelper(r);
			
		};
		
		
		/**
		 * @brief           Return raw data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void 
		GetRHelper              (Matrix< std::complex<float> >& M) {
			
			raw_data* rp = m_rrsi->rhelper();

			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = rp->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetSize(rp->dims); i++)
				M[i] = std::complex<float>(rp->dreal[i],rp->dimag[i]);
			
		};
		
		
		/**
		 * @brief           Set helper my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void
		SetHelper              (Matrix< double >& M) {
			
			helper_data h;

			h.dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				h.dims[j] = M.Dim(j);
			
			h.vals.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++)
				h.vals[i] = M[i];
			
			m_rrsi->helper(h);
			
		};
		
		
		/**
		 * @brief           Return helper data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void
		GetHelper              (Matrix< double >& M) {
			
			helper_data* hp = m_rrsi->helper();

			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = hp->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetSize(hp->dims); i++)
				M[i] = hp->vals[i];
			
		};
		
		
		/**
		 * @brief           Set kspace my data with ...
		 *
		 * @param  M        Given matrix
		 */
		void
		SetKSpace              (Matrix< double >& M) {
			
			helper_data h;

			h.dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				h.dims[j] = M.Dim(j);
			
			h.vals.length(M.Size());
			
			for (int i = 0; i < M.Size(); i++)
				h.vals[i] = M[i];
			
			m_rrsi->kspace(h);
			
		};
		
		
		/**
		 * @brief           Return kspace data after recon to ...
		 *
		 * @param  M        Given matrix
		 */
		void
		GetKSpace              (Matrix< double >& M) {
			
			helper_data* hp = m_rrsi->kspace();

			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = hp->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetSize(hp->dims); i++)
				M[i] = hp->vals[i];
			
		};
		
		
		/**
		 * @brief           Set my Pixel data
		 * 
		 * @param  M        Given matrix
		 */
		void 
		SetPixel            (Matrix<short>& M) {
			
			pixel_data p;

			p.dims.length(INVALID_DIM);
			
			for (int j = 0; j < INVALID_DIM; j++)
				p.dims[j] = M.Dim(j);
			
			p.vals.length(M.Size()); 
			
			for (int i = 0; i < M.Size(); i++)
				p.vals[i] = M[i];
			
			m_rrsi->pixel(p);
			
		};
		
		
		/**
		 * @brief           Put my Pixel data into 
		 * 
		 * @param  M        Given matrix
		 */
		void 
		GetPixel            (Matrix<short>& M) {
			
			pixel_data* pp = m_rrsi->pixel();

			for (int j = 0; j < INVALID_DIM; j++)
				M.Dim(j) = pp->dims[j];
			
			M.Reset();
			
			for (int i = 0; i < GetSize(pp->dims); i++)
				M[i] = pp->vals[i];
			
		};
		
		
		/**
		 * @brief           Raw repository size
		 *
		 * @return          Size
		 */
		long
		GetSize            (longs dims);
		
		
	private:
		
		RRSInterface_var    m_rrsi;       /**< Remote Recon interface               */
		CORBA::ORB_var      m_orb;        /**< Orb                                  */
		std::vector<short>  m_rstrats;    /**< Remote reconstruction strategies    */
		
	};
	
}


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
