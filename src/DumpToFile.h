#ifndef __DUMP_TO_FILE_H__
#define __DUMP_TO_FILE_H__

#include "ReconStrategy.h"

/**
 * @brief Dump data to file
 */
class DumpToFile : public ReconStrategy {


public:
	
	/**
	 * @brief Default constructor
	 */
	DumpToFile  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~DumpToFile () {};
	
	/**
	 * @brief Dump data to disk
	 */
	virtual RRSModule::error_code
	ProcessData ();
	
};

#endif /* __DUMMY_RECON_H__ */
