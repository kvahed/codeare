#ifndef __FFT_HPP__
#define __FFT_HPP__

#include "Matrix.hpp"

class FFT {
	

public:

	/**
	 * @brief   Forward transform
	 *
	 * @param   To transform
	 * @return  Transform
	 */
	static Matrix<cplx> 
	Forward     (const Matrix<cplx>& m);
	

	/**
	 * @brief   Backward transform
	 *
	 * @param   To transform
	 * @return  Transform
	 */
	static Matrix<cplx> 
	Backward    (const Matrix<cplx>& m);
	

	/**
	 * @brief   FFT shift
	 *
	 * @param   To shift
	 * @return  Shifted
	 */
	static Matrix<cplx> 
	Shift    (const Matrix<cplx>& m);
	

private:
	
	FFT()  {};
	~FFT() {};
	
};

#endif
