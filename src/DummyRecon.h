#ifndef __DUMMY_RECON_H__
#define __DUMMY_RECON_H__

#include "ReconStrategy.h"

/**
 * @brief Empty recon for test purposes
 */
class DummyRecon : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	DummyRecon  () {
		
	};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~DummyRecon () {};
	
	/**
	 * @brief Do nothing 
	 */
	virtual RRSModule::error_code
	ProcessData () { 

		/*for (int i = 0; i < m_raw.Size(); i++)
			m_raw[i] = m_raw[i];
		
		for (int i = 0; i < m_pixel.Size(); i++)
			m_pixel[i] = m_pixel[i];
		
		for (int i = 0; i < m_helper.Size(); i++)
		m_helper[i] = m_helper[i];*/
		
		return RRSModule::OK;

	};
	
};

#endif /* __DUMMY_RECON_H__ */
