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

#ifndef __NUFFT_H__
#define __NUFFT_H__

#include "ReconStrategy.h"

/**
 * @brief Non uniform FFT
 */
class NuFFT : public ReconStrategy {

public:
	
	/**
	 * @brief Default constructor
	 */
	NuFFT  ();
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~NuFFT () {};
	
	/**
	 * @brief Dump data to disk
	 */
	virtual RRSModule::error_code
	Process ();
	
private:
	
	int N [2]      = {0,0};
	int Nk[2]      = {0,0};
	
	// We only support 2 dims for the time being
	const int m_dim  = 2;
	
	// Some more stuff
	int       m_M    = 2;
	float     m_b    = 2.4;
	int       m_q    = 32;

	float*    m_tra;
	
	double* k_data;
	double* r_data;
	double* k_d;
	double* r_d;
	
};

#endif /* __NUFFT_H__ */
