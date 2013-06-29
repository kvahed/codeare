#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "BLACS.hpp"
#include "config.h"

#include <ostream>

/**
 * @brief BLACS grid 
 */
class Grid {
    
public:
    
	int  np; /**< @brief # of processes */
	int  rk; /**< @brief my rank        */
	int  ct; /**< @brief context        */
	int  nr; /**< @brief # of proc rows */
	int  nc; /**< @brief # of proc columns */
	int  mr; /**< @brief my row # */
	int  mc; /**< @brief my col # */
	char order; /**< @brief row/col major */
    
    const char* c_str() const;
    
    static Grid& Instance();

    virtual ~Grid();
    
private:
    
    Grid ();
    static Grid* m_inst;
    
};


/**
 * @brief            Dump to ostream
 *
 * @param  os        Output stream
 * @param  w         Workspace
 * @return           The output stream
 */
inline static std::ostream&
operator<< (std::ostream& os, Grid& g) {
	os << g.c_str();
    return os;
}

#endif //__GRID_HPP__
