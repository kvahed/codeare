#ifndef __DWT_HPP__
#define __DWT_HPP__

#include "Matrix.hpp"

class DWT {
	

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
	

private:
	
	DWT()  {};
	~DWT() {};


	/**
	 * @brief   Transform
	 *
	 * @param       To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	static Matrix<cplx> 
	Transform    (const Matrix<cplx>& m, bool bw);

	
};

#endif
