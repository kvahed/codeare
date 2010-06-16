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
	MedianFilter  () {m_ww = 25; m_wh = 25; };
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~MedianFilter () {};

	/**
	 * @brief Apply Median filter to image space
	 */
	virtual RRSModule::error_code ProcessData ();

	void SetFilterWidth (int val) { m_ww = val; };	
	void SetFilterHeight(int val) { m_wh = val; };	
	void SetFilterSize  (int val) { m_ww = val;  m_wh = val; };


private:	

	int m_ww; /**> @brief filter window width*/
	int m_wh; /**> @brief filter window height*/

};

#endif /* __MEDIAN_FILTER_H__ */
