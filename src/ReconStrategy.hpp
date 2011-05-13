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

#ifndef __RECON_STRATEGY_HPP__
#define __RECON_STRATEGY_HPP__

#include "Matrix.hpp"
#include "Configurable.hpp"


#include "DllExport.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif


#include <cstdlib>
#include <complex>


using namespace RRSModule;


namespace RRServer {
/**
 * @brief Strategy for reconstruction strategies
 *        Derive hereof to expand the reconstruction toolbox
 *
 */
class ReconStrategy : 
	public Configurable {


public:

	/**
	 * @brief Missing
	 */ 
	ReconStrategy  () {
		
	};
	
	/**
	 * @brief Missing
	 */ 
	virtual
	~ReconStrategy () {

	};


	/**
	 * @brief Data procession function 
	 */ 
	virtual error_code
	Process     () = 0;
	
	/**
	 * @brief Initilise
	 */ 
	virtual error_code
	Init     () {};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetRaw           (raw_data* raw)   {
		
		for (int j = 0; j < INVALID_DIM; j++)
			raw->dims[j] = m_raw.Dim(j);
			
		raw->dreal.length(m_raw.Size()); 
		raw->dimag.length(m_raw.Size());
			
		for (int i = 0; i < m_raw.Size(); i++) {
			raw->dreal[i] = m_raw[i].real();
			raw->dimag[i] = m_raw[i].imag(); 
		}
		
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRaw           (const raw_data* raw)   {

		for (int i = 0; i < INVALID_DIM; i++)
			m_raw.Dim(i) = raw->dims[i];
		
		m_raw.Reset ();
		
		for (int j = 0; j < m_raw.Size(); j++)
			m_raw[j] =  std::complex<float> (raw->dreal[j], raw->dimag[j]);
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetRHelper           (raw_data* rhelper)   {
		
		for (int j = 0; j < INVALID_DIM; j++)
			rhelper->dims[j] = m_rhelper.Dim(j);
		
		rhelper->dreal.length(m_rhelper.Size()); 
		rhelper->dimag.length(m_rhelper.Size());
		
		for (int i = 0; i < m_rhelper.Size(); i++) {
			rhelper->dreal[i] = m_rhelper[i].real();
			rhelper->dimag[i] = m_rhelper[i].imag(); 
		}
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRHelper           (const raw_data* rhelper)   {

		for (int i = 0; i < INVALID_DIM; i++)
			m_rhelper.Dim(i) = rhelper->dims[i];
		
		m_rhelper.Reset ();
		
		for (int j = 0; j < m_rhelper.Size(); j++)
			m_rhelper[j] =  std::complex<float> (rhelper->dreal[j], rhelper->dimag[j]);
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetHelper           (helper_data* helper)   {

		for (int i = 0; i < m_helper.Size(); i++)
			helper->vals[i] = m_helper[i];

	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetHelper           (const helper_data* helper)   {

		for (int i = 0; i < INVALID_DIM; i++)
			m_helper.Dim(i) = helper->dims[i];

		m_helper.Reset ();

		for (int j = 0; j < m_helper.Size(); j++)
			m_helper[j] =  helper->vals[j];
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetKSpace           (helper_data* kspace)   {

		for (int i = 0; i < m_kspace.Size(); i++)
			kspace->vals[i] = m_kspace[i];

	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetKSpace           (const helper_data* kspace)   {

		for (int i = 0; i < INVALID_DIM; i++)
			m_kspace.Dim(i) = kspace->dims[i];

		m_kspace.Reset ();

		for (int j = 0; j < m_kspace.Size(); j++)
			m_kspace[j] =  kspace->vals[j];
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetPixel           (pixel_data* pixel)   {

		for (int i = 0; i < m_pixel.Size(); i++)
			pixel->vals[i] = m_pixel[i];

	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetPixel           (const pixel_data* pixel)   {

		for (int i = 0; i < INVALID_DIM; i++)
			m_pixel.Dim(i) = pixel->dims[i];

		m_pixel.Reset ();

		for (int j = 0; j < m_pixel.Size(); j++) 
			m_pixel[j] =  pixel->vals[j];
		
	};
	
protected:

	Matrix<raw>     m_raw;         /*!< raw data matrix                    */
	Matrix<raw>     m_rhelper;     /*!< raw helper matrix                  */
	Matrix<double>  m_helper;      /*!< helper matrix                      */
	Matrix<double>  m_kspace;      /*!< kspace matrix                      */
	Matrix<short>   m_pixel;       /*!< pixel data matrix                  */

};



#endif /* __RECON_STRATEGY_HPP__ */

}

typedef RRServer::ReconStrategy* create_t  ();
typedef void           destroy_t (RRServer::ReconStrategy*);

