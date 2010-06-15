#ifndef __RECON_SERVANT_H__
#define __RECON_SERVANT_H__

#include "Matrix.h"
#include "RRSModule.hh"

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
