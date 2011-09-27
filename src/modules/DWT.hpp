#ifndef __DWT_HPP__
#define __DWT_HPP__

#include "Matrix.hpp"

/**
 * @brief 2D Discrete wavelet transform for Matrix template<br/>(Daubechies wavelets)
 */
class DWT {
	

public:

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	static Matrix<cplx> 
	Forward     (const Matrix<cplx>& m);
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	static Matrix<cplx> 
	Backward    (const Matrix<cplx>& m);
	

private:
	
	/**
	 * Static class
	 */
	DWT()  {};

	/**
	 * Static class
	 */
	~DWT() {};

	/**
	 * @brief   Transform
	 *
	 * @param   m   To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	static Matrix<cplx> 
	Transform    (const Matrix<cplx>& m, bool bw);

	
};

#endif
