/*
 * NLCG.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _NLCG_HPP_
#define _NLCG_HPP_

#include <NonLinear.hpp>

#ifdef _MSC_VER
std::string ofstr = "    %02Iu - nrms: %1.4e, l-search: %d, ";
#else
std::string ofstr = "    %02zu - nrms: %1.4e, l-search: %d, ";
#endif

namespace codeare {
    namespace optimisation {
        
template<class T> class NLCG : public NonLinear<T> {

    typedef typename TypeTraits<T>::RT real_t;

public:
    
    NLCG (const Params& p) : NonLinear<T>::NonLinear (p) {
        _nliter = try_to_fetch<int> (p, "nliter", 0);
        _cgconv = try_to_fetch<float> (p, "cgconv", 0);
        _lsiter = try_to_fetch<int> (p, "lsiter", 0);
        _lsa = try_to_fetch<float> (p, "lsa", 0);
        _lsb = try_to_fetch<float> (p, "lsb", 0);
        _pls = try_to_fetch<bool> (p, "parallel_linesearch", false);
    }
    
	virtual ~NLCG() {};

    inline virtual size_t LineSearch (const Operator<T>* A, const Matrix<T>& x, const real_t& t0,
                                      const real_t& f0, real_t& rmse, real_t& t) const {
        size_t i = 0;
        real_t f1 = 0, g0dx = abs(_g0.dotc(_dx));
        while (i < _lsiter) {
            t = t0 * pow(_lsb,i+1);
            f1 = A->obj (x, _dx, t, rmse);
            if (f1 <= f0 - _lsa * t * g0dx)
                break;
            ++i;
        }
        return i;
    }

    inline virtual size_t LineSearchParallel (const Operator<T>* A, const Matrix<T>& x, const real_t& t0,
                                              const real_t& f0, real_t& rmse, real_t& t) const {
        Vector<real_t> rmses (_lsiter);
        Vector<real_t> ts (_lsiter);
        Vector<real_t> f1s (_lsiter);
        real_t g0dx = abs(_g0.dotc(_dx));
        int li = -1;
#pragma omp parallel for default(shared) schedule(static,1) num_threads(_lsiter)
        for (size_t i = 0; i < _lsiter; ++i) {
            ts[i] = t0 * pow(_lsb,i+1);
            f1s[i] = A->obj (x, _dx, ts[i], rmses[i]);
        }
        for (size_t i = 0; i < _lsiter; ++i)
            if (f1s[i] <= f0 - _lsa * t * g0dx) {
                li = i;
                break;
            }
        t = ts[li];
        rmse = rmses[li];
        return (li>-1) ? li : _lsiter;
    }

    inline virtual void Minimise (Operator<T>* A, Matrix<T>& x) {

        real_t t0  = 1.0, t = 1.0, z = 0., xn = norm(x), rmse, bk, f0, dxn;
        Vector<real_t> rms(_lsiter);
        Vector<size_t> pos(_lsiter);

        _g0  = A->df (x);
        _dx  = -_g0;
    
        for (size_t k = 0; k < _nliter; k++) {
        
            A->Update(_dx);
        
            f0 = A->obj (x, _dx, z, rmse);

            size_t li = _pls ?
                LineSearchParallel (A, x, t0, f0, rmse, t) :
                LineSearch (A, x, t0, f0, rmse, t);
            printf (ofstr.c_str(), k, rmse, li); fflush (stdout);
            if (li == _lsiter) {
                printf ("Reached max line search, exiting... \n");
                return;
            }
        
            if      (li > 2) t0 *= _lsb;
            else if (li < 1) t0 /= _lsb;
        
            // Update image
            x  += _dx * t;
        
            // CG computation
            _g1 =  A->df (x);
            bk  =  real(_g1.dotc(_g1)) / real(_g0.dotc(_g0));
            _g0 =  _g1;
            _dx = -_g1 + _dx * bk;
            dxn =  norm(_dx)/xn;
        
            printf ("dxnrm: %0.4f\n", dxn);
            if (dxn < _cgconv)
                break;
        
        }

    }


    virtual std::ostream& Print (std::ostream& os) const {
		NonLinear<T>::Print(os);
        os << "    Iterations: (" << _nliter << ") LS(" << _lsiter << ")" << std::endl;
        os << "    Conv: CG(" << _cgconv << ")" << std::endl;
        os << "    LS (" << (_pls ? "parallel" : "linear") << ")" << std::endl; 
        os << "    LS brackets: lsa(" << _lsa << ") lsb(" << _lsb << ")";
		return os;
    }
    
private:    

    Matrix<T> _dx, _g0, _g1;
    bool _pls;
    typename TypeTraits<T>::RT _cgconv, _lsa, _lsb;
    size_t _lsiter, _nliter;
    

};
        
    }}

#endif // _NLCG_HPP_
