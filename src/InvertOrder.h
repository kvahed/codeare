#ifndef __INVERT_ORDER_H__
#define __INVERT_ORDER_H__

#include "ReconStrategy.h"

/**
 * @brief Reorder for test purposes
 */
class InvertOrder : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	InvertOrder  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~InvertOrder () {};
	
	/**
	 * @brief Invert the order of the data. 
	 *        No good for anything besides testing.
	 */
	virtual RRSModule::error_code
	ProcessData ();
	
};

#endif /* __INVERT_ORDER_H__ */
