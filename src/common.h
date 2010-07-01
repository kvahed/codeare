#ifndef __COMMON_H__
#define __COMMON_H__

/**
 * Method dispatch enumeration
 */
enum      method       {
	
	DUMMY_RECON,
	HANN_WINDOW,
	MEDIAN_FILTER,
	INVERT_ORDER,
	DUMP_TO_FILE,
	MEDIAN_FILTER_OMP
	
};


/**
 * Returned error codes 
 */
enum      error_code   {
	
	OK,
	ERROR_GENERAL,
	ERROR_UNKNOWN_METHOD
	
};


#endif //__COMMON_H__
