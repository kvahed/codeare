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

#ifndef __RECON_SERVANT_H__
#define __RECON_SERVANT_H__

#include "Matrix.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif


using namespace RRSModule;


/**
 * @brief Servant implementation 
 *        Perform duties on remote server
 */
class ReconServant : 
	public POA_RRSModule::RRSInterface , 
	public PortableServer::RefCountServantBase {
	
public:
	
	/**
	 * @brief Default constructor
	 */
	ReconServant  ();
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~ReconServant ();

	/**
	 * @brief Dispatch data to available methods 
	 */
	virtual error_code
	process_data  (method m);
	
	virtual error_code 
	Process (const char* name);

	/**
	 * @brief Get data from recon
	 */
	raw_data* 
	raw           ();
	
	/**
	 * @brief Set data for recon
	 */
	void 
	raw           (const raw_data&);
	
	/**
	 * @brief Get data from recon
	 */
	raw_data* 
	helper        ();
	
	/**
	 * @brief Set data for recon
	 */
	void 
	helper        (const raw_data&);
	
	/**
	 * @brief Get data from recon
	 */
	pixel_data* 
	pixel         ();
	
	/**
	 * @brief Set data for recon
	 */
	void 
	pixel         (const pixel_data&);
	
	/**
	 * @brief Get data from recon
	 */
	strings* 
	labels        ();
	
	/**
	 * @brief Set data for recon
	 */
	void 
	labels        (const strings&);
	
	
private:
	
	raw_data      m_raw;    /**< Raw data    (complex float) */
	raw_data      m_helper; /**< Image data  (complex float) */
	pixel_data    m_pixel;  /**< Helper data (short)         */
	strings       m_labels;

	bool          m_have_raw;
	bool          m_have_helper;
	bool          m_have_pixel;
	bool          m_have_labels;

};

#endif // __RECON_SERVANT_H__
