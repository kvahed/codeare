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

#ifndef __RECON_STRATEGY_H__
#define __RECON_STRATEGY_H__

#include "Matrix.h"
#include "Configurable.h"


#include "DllExport.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif


#include <cstdlib>
#include <complex>


using namespace RRSModule;

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
	 * @brief Get data from recon
	 */
	void 
	GetRaw           (raw_data* raw)   {
		
		for (int i = 0; i < m_raw.Size(); i++) {
			raw->dreal[i] = abs(m_raw[i]);
			raw->dimag[i] = arg(m_raw[i]); 
		}
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRaw           (const raw_data* raw)   {

		int dim[16];
		int i = 0;

		for (i = 0; i < 16; i++)
			dim[i] = raw->dims[i];

		m_raw.Reset (dim);

		for (i = 0; i < m_raw.Size(); i++)
			m_raw[i] =  std::complex<float> (raw->dreal[i], raw->dimag[i]);
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetHelper           (raw_data* helper)   {

		for (int i = 0; i < m_helper.Size(); i++) {
			helper->dreal[i] = abs(m_helper[i]);
			helper->dimag[i] = arg(m_helper[i]); 
		}
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetHelper           (const raw_data* helper)   {

		int dim[16];
		int i = 0;

		for (i = 0; i < 16; i++)
			dim[i] = helper->dims[i];

		m_helper.Reset (dim);

		for (i = 0; i < m_helper.Size(); i++)
			m_helper[i] =  std::complex<float> (helper->dreal[i], helper->dimag[i]);
		
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

		int dim[16];
		int i = 0;

		for (i = 0; i < 16; i++)
			dim[i] = pixel->dims[i];

		m_pixel.Reset (dim);

		for (i = 0; i < m_pixel.Size(); i++) 
			m_pixel[i] =  pixel->vals[i];
		
	};
	
protected:

	Matrix< std::complex<float> > m_raw;         /*!< raw data matrix                    */
	Matrix< std::complex<float> > m_helper;      /*!< helper matrix                      */
	Matrix< short >               m_pixel;       /*!< pixel data matrix                  */

};

#endif /* __RECON_STRATEGY_H__ */

typedef ReconStrategy* create_t  ();
typedef void           destroy_t (ReconStrategy*);
