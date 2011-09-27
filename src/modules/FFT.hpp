#ifndef __FFT_HPP__
#define __FFT_HPP__

#include "Matrix.hpp"

/**
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */
class FFT {

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
	

	/**
	 * @brief    FFT shift
	 *
	 * @param  m To shift
	 * @return   Shifted
	 */
	static Matrix<cplx> 
	Shift    (const Matrix<cplx>& m);
	

private:
	
	/**
	 * Static class
	 */
	FFT()  {};


	/**
	 * Static class
	 */
	~FFT() {};
	
};

#endif
