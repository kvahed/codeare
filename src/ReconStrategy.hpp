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
#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

#include "cycle.h"            // FFTW cycle implementation

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
	 * @brief Initilise
	 */ 
	virtual error_code
	Finalise     () {};
	
	/**
	 * @brief  CPU clock rate
	 *
	 * @return clock rate
	 */
	
	double 
	ClockRate () {

#if defined(HAVE_MACH_ABSOLUTE_TIME)

		// OSX

		uint64_t freq = 0;
		size_t   size = sizeof(freq);
		
		if (sysctlbyname("hw.tbfrequency", &freq, &size, NULL, 0) < 0)
			perror("sysctl");
		
		return freq;

#else
		
		// LINUX

		FILE* scf;
		std::string fname = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq";
		int   freq = 1;
		
		scf = fopen(fname.c_str(), "rb");
		
		if (scf != NULL) {
			int read = fscanf(scf,"%i",&freq);
#ifdef VERBOSE
			printf ("Read %i from %s\n", read, fname.c_str());
#endif
			fclose(scf);
		}
		
		return 1000.0 * freq;
		
#endif
	}


	/**
	 * @brief Get data from recon
	 */
	void
	GetRaw           (raw_data* r)   {
		
		for (int j = 0; j < INVALID_DIM; j++)
			r->dims[j] = m_raw.Dim(j);
			
		r->dreal.length(m_raw.Size()); 
		r->dimag.length(m_raw.Size());
			
		for (int i = 0; i < m_raw.Size(); i++) {
			r->dreal[i] = m_raw[i].real();
			r->dimag[i] = m_raw[i].imag(); 
		}

		m_raw.Clear();

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRaw           (const raw_data* r)   {
		
		for (int i = 0; i < INVALID_DIM; i++)
			m_raw.Dim(i) = r->dims[i];
		
		m_raw.Reset ();
		
		for (int j = 0; j < m_raw.Size(); j++)
			m_raw[j] =  std::complex<float> (r->dreal[j], r->dimag[j]);
		
	}
	
	/**
	 * @brief Get data from recon
	 */
	void
	GetRaw           (Matrix<raw>* m)   {

		(*m) = m_raw;
		m_raw.Clear();

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRaw           (const Matrix<raw>* m)   {
		
		m_raw = (*m);
		
	}
	
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
	GetRHelper           (Matrix<raw>* m)   {

		m = &m_rhelper;

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRHelper           (const Matrix<raw>* m)   {
		
		m_rhelper = (*m);
		
	}
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetHelper           (helper_data* helper)   {

		for (int j = 0; j < INVALID_DIM; j++)
			helper->dims[j] = m_helper.Dim(j);
		
		helper->vals.length(m_helper.Size()); 
		
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
	GetHelper           (Matrix<double>* m)   {

		m = &m_helper;

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetHelper           (const Matrix<double>* m)   {
		
		m_helper = (*m);
		
	}
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetKSpace           (helper_data* kspace)   {

		for (int j = 0; j < INVALID_DIM; j++)
			kspace->dims[j] = m_kspace.Dim(j);
		
		kspace->vals.length(m_kspace.Size()); 
		
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
		
	}
	
	/**
	 * @brief Get data from recon
	 */
	void
	GetKSpace           (Matrix<double>* m)   {

		m = &m_kspace;

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetKSpace           (const Matrix<double>* m)   {
		
		m_kspace = (*m);
		
	}
	
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

	/**
	 * @brief Get data from recon
	 */
	void
	GetPixel           (Matrix<short>* m)   {

		m = &m_pixel;

	}

	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetPixel           (const Matrix<short>* m)   {
		
		m_pixel = (*m);
		
	}
	


	void Name (const char* name) { m_name = std::string(name);}

	const char* Name () {return m_name.c_str();}
	
protected:

	Matrix<raw>     m_raw;         /*!< raw data matrix                    */
	Matrix<raw>     m_rhelper;     /*!< raw helper matrix                  */
	Matrix<double>  m_helper;      /*!< helper matrix                      */
	Matrix<double>  m_kspace;      /*!< kspace matrix                      */
	Matrix<short>   m_pixel;       /*!< pixel data matrix                  */

	std::string     m_name;

};



#endif /* __RECON_STRATEGY_HPP__ */

}

typedef RRServer::ReconStrategy* create_t  ();
typedef void           destroy_t (RRServer::ReconStrategy*);

