#ifndef __GENERIC_REGRID_H__
#define __GENERIC_REGRID_H__

#include "ReconStrategy.h"

/**
 * @brief Generic regridder
 *        Expects to identically sized matrices (data & k-space positions)
 */
class GenericRegrid : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	GenericRegrid  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~GenericRegrid () {};
	
	/**
	 * @brief Regrid data to Cartesian k-space
	 */
	virtual RRSModule::error_code
	ProcessData () {
		return RRSModule::OK;
	};
	
};

#endif /* __GENERIC_REGRID_H__ */
