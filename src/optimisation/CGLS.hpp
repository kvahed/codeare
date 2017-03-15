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

template<class T> class CGLS : public Linear<T> {

    typedef typename TypeTraits<T>::RT RT;

public:
    CGLS (const size_t& maxit = 10, const RT& epsilon = 1.0e-6,
          const RT& lambda = 1.0e-6, const int& verbosity = 0) :
        Linear<T>::Linear(verbosity), _verbosity(verbosity), _nrows(1),
        _ncols(1), _maxit(maxit), _epsilon(epsilon), _lambda(lambda) {}
    
    virtual ~CGLS () {}

    inline virtual Matrix<T> Solve (const Operator<T>& A, const MatrixType<T>& x) {
        _p = A/x;
        if (_maxit == 0)
            return _p;
        _r  = _p;
    _rn = _xn = std::real(dotc(_r,_p));
    Matrix<T> ret;
    for (size_t i = 0; i < _maxit; i++) {
      _res.push_back(_rn/_xn);
      if (i==0)
                ret = zeros<T>(size(_p));
            if (boost::math::isnan(_res[i]) || _res[i] <= _epsilon) {
                printf ("    %03zu %.7f\n", i, _res[i]);
                break;
      }
            if (_verbosity)
                printf ("    %03zu %.7f\n", i, _res[i]);
      _q  = A/(A*_p);
            if (_lambda)
                _q += _lambda * _p;
            _ts  = _rn / std::real(dotc(_p,_q));
            ret += _ts * _p;
            _r  -= _ts * _q;
            _rno = _rn;
            _rn  = std::real(dotc(_r,_r));
            _p  *= _rn / _rno;
            _p  += _r;
        }
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
