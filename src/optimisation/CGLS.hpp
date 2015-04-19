/*
 * CGLS.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */

#ifndef SRC_OPTIMISATION_CGLS_HPP_
#define SRC_OPTIMISATION_CGLS_HPP_

#include <Linear.hpp>
#include <Lapack.hpp>
#include <CX.hpp>

namespace codeare {
namespace optimisation {

template<class T>
class CGLS : public Linear<T> {
	typedef typename TypeTraits<T>::RT RT;
public:
	CGLS (const size_t& maxit = 10, const RT& epsilon = 1.0e-6,
			const RT& lambda = 1.0e-6, const int& verbosity = 0) :
				Linear<T>::Linear(verbosity), _verbosity(verbosity), _nrows(1),
				_ncols(1), _maxit(maxit), _epsilon(epsilon), _lambda(lambda) {}

	virtual ~CGLS () {}

	inline Matrix<T> Solve (const Operator<T>& A, const Matrix<T>& x) {

		Matrix<T> ret;
        Vector<Matrix<T> > vc;

        typedef typename Vector<T>::iterator it_type;

		_p = (A/x);// * m_ic;
		if (_maxit == 0)
			return _p;

		_r  = _p;
        _xn = std::real(_p.dotc(_p));
        _rn = _xn;

		for (size_t i = 0; i < _maxit; i++) {
			_res.push_back(_rn/_xn);
			if (i==0)
				ret = zeros<T>(size(_p));
			if (boost::math::isnan(_res[i]) || _res[i] <= _epsilon)
				break;
			if (_verbosity)
				printf ("    %03lu %.7f\n", i, _res[i]);
			_q  = A/(A*_p);
			if (_lambda)
				_q += _lambda * _p;
			_ts  = _rn / std::real(_p.dotc(_q));
			ret += _ts * _p;
			_r  -= _ts * _q;
			_rno = _rn;
			_rn  = std::real(_r.dotc(_r));
			_p  *= _rn / _rno;
			_p  += _r;
			if (_verbosity)
				vc.push_back(ret);
		}

		/*
        if (m_verbose) { // Keep intermediate results
            size_t cpsz = numel(ret);
            ret = Matrix<T> (size(ret,0), size(ret,1), (m_nx[0] == 3) ?
            		size(x,2) : 1, vc.size());
            it_type it = x.Begin();
            for (size_t i = 0; i < vc.size(); i++) {
                std::copy (vc[i].Begin(), vc[i].End(), it);
                it += cpsz;
            }
            vc.Clear();
        }
*/
		return ret;// * m_ic;
	}

protected:
	size_t _nrows, _ncols, _maxit;
	RT _epsilon, _rel_mat_err, _rel_rhs_err, _lambda, _ts;
	RT _rn, _xn, _rno;
	Matrix<T> _p, _r, _q;
    Vector<Matrix<cxfl> > vc;
    Vector<RT> _res;
	int _verbosity;

};
}}

#endif /* SRC_OPTIMISATION_CGLS_HPP_ */
