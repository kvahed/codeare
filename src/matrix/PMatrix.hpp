#ifndef __P_MATRIX__
#define __P_MATRIX__

#include "ScalapackTraits.hpp"
#include "Matrix.hpp"

static int izero = 0;

/**
 * @brief MPI aware C++ friendly matrix.<br/>
 *        This is now only 2D for Scalapack operations only.<br/>
 *        Maybe there is more general need for higher dimensions later.
 * 
 */
template <class T> 
class PMatrix : public Matrix<T> {



 public: 

	/**
	 * @brief             Default constructor
	 */
	PMatrix () : _bs (16) {

		// Validate
		T t;
		PMatrix<T>::Validate (t);

		// Get at home
		GridInfo();
		
	}
	

	/**
	 * @brief             Construct with sizes
	 */
	PMatrix (const size_t& cols, const size_t& rows) : _bs (16) {

		// Validate
		T t;
		PMatrix<T>::Validate (t);

		int info; 

		// Get at home
		GridInfo();

		// Global size
		_gdim[0] = cols;
		_gdim[1] = rows;
		
		// Local size (only with MPI different from global)
#ifdef HAVE_MPI
		Matrix<T>::_dim[0] = numroc_ (&_gdim[0], &_bs, &_gd.mc, &izero, &_gd.nc);
		Matrix<T>::_dim[1] = numroc_ (&_gdim[1], &_bs, &_gd.mr, &izero, &_gd.nr);
#else
		Matrix<T>::_dim[0] = _gdim[0];
		Matrix<T>::_dim[1] = _gdim[1];
#endif
		
		// Allocate
		Matrix<T>::_M.resize(Matrix<T>::Size());
		
		// Descriptor 
#ifdef HAVE_MPI
		int dims[2]; dims[0] = Matrix<T>::_dim[0]; dims[1] = Matrix<T>::_dim[1];
		descinit_(_desc, &_gdim[0], &_gdim[1], &_bs, &_bs, &izero, &izero, &_gd.ct, 
				  dims, &info);
#endif
		
#ifdef DESC_DEBUG
		printf ("info(%d) desc({%d, %d, %4d, %4d, %d, %d, %d, %d, %4d})\n", 
				info,     _desc[0], _desc[1], _desc[2], _desc[3], 
				_desc[4], _desc[5], _desc[6], _desc[7], _desc[8]);
#endif


	}
	

	/**
	 * @brief          Virtual destructor
	 */
	virtual 
	~PMatrix () {}


	/**
	 * @brief         Get # of elements
	 *
	 * @return        # of elements
	 */
	inline size_t 
	GHeight () const {
		return _gdim[0];
	}
	
	
	/**
	 * @brief         Get # of elements
	 *
	 * @return        # of elements
	 */
	inline size_t 
	GWidth () const {
		return _gdim[1];
	}
	
	
	/**
	 * @brief         Get # of elements
	 *
	 * @return        # of elements
	 */
	inline size_t 
	g_m () const {
		return _gdim[0];
	}
	
	
	/**
	 * @brief         Get # of elements
	 *
	 * @return        # of elements
	 */
	inline size_t 
	g_n () const {
		return _gdim[1];
	}
	
	
	/**
	 * @brief         Get descriptor vector
	 *
	 * @return        Descriptor vector
	 */
	inline const int*
	Desc () const {
		return _desc;
	}
	

 protected:

	/**
	 * @brief   Who are we and where are we?
	 */
	inline void 
	GridInfo () {

		// Defaults
		_gd.rk = 0; 
		_gd.np = 0; 
		_gd.mr = 0; 
		_gd.mc = 0; 
		_gd.nc = 0; 
		_gd.nr = 0; 
		_gd.ct = 0;

#ifdef HAVE_MPI		
		Cblacs_pinfo    (&_gd.rk, &_gd.np);
		Cblacs_get      (-1, 0, &_gd.ct);
		Cblacs_gridinfo (_gd.ct, &_gd.nr, &_gd.nc, &_gd.mr, &_gd.mc); 
#endif

#ifdef GINFO_DEBUG
		printf("id (%d), row(%d), col(%d)\n", _gd.rk, _gd.mr, _gd.mc);
#endif

	}


	// BLACS 
	grid_dims        _gd;
	int              _desc[9]; /**< @brief matrix grid vector */
	int              _bs;
	
	// Data
	int              _gdim[2]; /**< @brief Global dimensions */
	//int              _dim[2];
	//std::valarray<T> _M;

};

#endif //__P_MATRIX__
