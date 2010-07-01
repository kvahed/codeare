#ifndef __MEDIAN_FILTER_OMP_H__
#define __MEDIAN_FILTER_OMP_H__

#include "ReconStrategy.h"

/**
 * @brief Median filter with OpenMP support
 */
class MedianFilter_OMP : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	MedianFilter_OMP  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~MedianFilter_OMP () {};

	/**
	 * @brief Apply Median filter to image space
	 */
	virtual RRSModule::error_code
	ProcessData ();

};

#endif /* __MEDIAN_FILTER_OMP_H__ */
