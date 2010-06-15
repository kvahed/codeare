#ifndef __HANN_WINDOW_H__
#define __HANN_WINDOW_H__

#include "ReconStrategy.h"

/**
 * @brief Hann window
 */
class HannWindow : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	HannWindow  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~HannWindow () {};
	
	/**
	 * @brief Apply Hann window to k-space
	 */
	virtual RRSModule::error_code
	ProcessData () {
		return RRSModule::OK;
	};
	
};

#endif /* __HANN_WINDOW_H__ */
