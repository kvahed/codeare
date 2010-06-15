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
	DummyRecon  () {};
	
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
		return RRSModule::OK;
	};
	
};

#endif /* __DUMMY_RECON_H__ */
