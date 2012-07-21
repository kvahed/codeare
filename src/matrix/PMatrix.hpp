#ifndef __P_MATRIX__
#define __P_MATRIX__

#include "ScalapackTraits.hpp"

#include <valarray>

static int izero = 0;

/**
 * @brief MPI aware C++ friendly matrix.<br/>
 *        This is now only 2D
 * 
 */
template <class T> 
class PMatrix {

 public: 

	PMatrix () : _bs (16) {
		
		// Get at home
		GridInfo();
		
	}
	
	PMatrix (const size_t& cols, const size_t& rows) : _bs (16) {

		// Get at home
		GridInfo();

		// Global size
		_gdim[0] = cols;
		_gdim[1] = rows;
		
		// Local size
		_dim[0] = numroc_ (&_gdim[0], &_bs, &_gd.mc, &izero, &_gd.nc);
		_dim[1] = numroc_ (&_gdim[1], &_bs, &_gd.mr, &izero, &_gd.nr);
		
		// Allocate
		_M.resize(Size());
		
	}
	
	inline size_t 
		Size() const {
		return _dim[0] * _dim[1];
	}
	
	
 private:

	/**
	 * @brief   Who are we and where are we?
	 */
	inline void GridInfo () {

		// Defaults
		_gd.rk = 0; 
		_gd.np = 0; 
		_gd.mr = 0; 
		_gd.mc = 0; 
		_gd.nc = 0; 
		_gd.nr = 0; 
		_gd.ct = 0;
		
		Cblacs_pinfo    (&_gd.rk, &_gd.np);
		Cblacs_get      (-1, 0, &_gd.ct);
		Cblacs_gridinfo (_gd.ct, &_gd.nr, &_gd.nc, &_gd.mr, &_gd.mc); 

#ifdef MPI_DEBUG
		printf("id (%d), row(%d), col(%d)\n", _gd.rk, _gd.mr, _gd.mc);
#endif

	}

	// BLACS 
	grid_dims        _gd;
	int              _desc[9]; /**< @brief matrix grid vector */
	int              _bs;
	
	// Data
	std::valarray<T> _M;      /**< @brief The meat */
	int              _dim[2]; /**< @brief Local dimensions */
	int              _gdim[2]; /**< @brief Global dimensions */

};

#endif //__P_MATRIX__
