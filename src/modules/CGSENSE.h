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

#ifndef __CGSENSE_H__
#define __CGSENSE_H__

#include "ReconStrategy.h"

#include <noncart/nufft.h>

typedef std::complex<float> raw;

/**
 * @brief Non uniform FFT
 */
class CGSENSE : public ReconStrategy {


public:
	
	/**
	 * @brief Default constructor
	 */
	CGSENSE  ();
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~CGSENSE () {};
	
	/**
	 * @brief Dump data to disk
	 */
	virtual RRSModule::error_code
	Process ();
	
 

private:

	/**
	 * @brief Multiply with Hermitian counterpart of the k-space sampling 
	 */
	void
	EH (Matrix< raw >* in, Matrix< raw >* out);

	/**
	 * @brief Multiply with spatial image data of 
	 */
	void
	E  (Matrix< raw >* in, Matrix< raw >* out);

	Matrix < raw >     m_sens;
	Matrix < raw >     m_temp;

	unsigned short int m_iter;
	
	int                m_verbose;

	noncart::nufft     m_nufft;

};

#endif /* __CGSENSE_H__ */
