#ifndef __MEDIAN_FILTER_H__
#define __MEDIAN_FILTER_H__

#include "ReconStrategy.h"

/**
 * @brief Median filter
 */
class MedianFilter : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	MedianFilter  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~MedianFilter () {};

	/**
	 * @brief Apply Median filter to image space
	 */
	virtual RRSModule::error_code
	ProcessData () {
		return RRSModule::OK;
	};
		
};

#endif /* __MEDIAN_FILTER_H__ */
