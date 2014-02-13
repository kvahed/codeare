/*
 * AssignmentHandler.hpp
 *
 *  Created on: Nov 21, 2013
 *      Author: kvahed
 */

#ifndef ASSIGNMENTHANDLER_HPP_
#define ASSIGNMENTHANDLER_HPP_

template<class T>
class AssignmentHandler {

public:
	inline AssignmentHandler (Matrix<T>& M) : _M(M) {

	}

	inline AssignmentHandler& operator= (const Matrix<T>& M) {
		_M = M;
		return *this;
	}

protected:
	Matrix<T>& _M;
};


#endif /* ASSIGNMENTHANDLER_HPP_ */
