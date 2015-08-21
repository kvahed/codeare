/*
 * LBFGS.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _LBFGS_HPP_
#define _LBFGS_HPP_

#include <NonLinear.hpp>

#define LBFGS_FLOAT 32
#include "lbfgs.h"


namespace codeare {
    namespace optimisation {
        
template<class T>
class LBFGS : public NonLinear<T> {

    typedef typename TypeTraits<T>::RT RT;
    
public:
    LBFGS () {
        lbfgs_parameter_init(&_p);
    }
	LBFGS (const Params& p) : NonLinear<T>::NonLinear(p)  {
        lbfgs_parameter_init(&_p);
        _p.max_iterations = try_to_fetch<int> (p, "nliter", 0);
        _p.min_step = try_to_fetch<float> (p, "lsa", 0);
        _p.max_step = try_to_fetch<float> (p, "lsb", 0);
    }
    LBFGS (const LBFGS& tocopy) {};
    virtual ~LBFGS () {}

    inline static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step, int n, int k, int ls) {
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    inline static lbfgsfloatval_t evaluate (void* instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
            const int n, const lbfgsfloatval_t step) {
        std::copy (x, x+n, (lbfgsfloatval_t*)&_x[0]);
        _g = _A->df (_x);
        (_A)->Update(_g);
        std::copy ((lbfgsfloatval_t*)&_g[0], (lbfgsfloatval_t*)&g[0]+n, &_x[0]);
        return (lbfgsfloatval_t) _A->obj(_x, _g, step, _rmse) / (lbfgsfloatval_t)n;
    }
    
    inline virtual void Minimise (Operator<T>* A, Matrix<T>& x) {
        _A = A;
        _x = x;
        _g = x;
        lbfgsfloatval_t fx;
        int ret = lbfgs(2*x.Size(), (lbfgsfloatval_t*)&x[0], &fx, evaluate, progress, NULL, &_p);
        printf("    L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

    }
    virtual std::ostream& Print(std::ostream& os) const {
        NonLinear<T>::Print(os);
        os << "    m(" << _p.m << ") epsilon(" << _p.epsilon << ") past(" << _p.past << ") delta("
           << _p.delta << ")" << std::endl;
        os << "    max_iterations(" << _p.max_iterations << ") linesearch(" << _p.linesearch
           << ") max_linesearch(" << _p.max_linesearch << ")" << std::endl;
        os << "    min_step(" << _p.min_step << ") max_step(" << _p.max_step << ")" << std::endl;
        os << "    ftol(" << _p.ftol << ") wolfe(" << _p.wolfe << ") gtol(" << _p.gtol
           << ") xtol(" << _p.xtol << ")" << std::endl;
        os << "    orthantwise_c(" << _p.orthantwise_c << ") orthantwise_start("
           << _p.orthantwise_start << ") orthantwise_end(" << _p.orthantwise_end << ")" << std::endl;
        return os;
    }

    lbfgs_parameter_t _p;
    static Matrix<T> _x;
    static Matrix<T> _g;
    static RT _rmse;
    static Operator<T>* _A;
};

        template<class T> Operator<T>* LBFGS<T>::_A = 0;
        template<class T> Matrix<T> LBFGS<T>::_x = Matrix<T>(256,256);
        template<class T> Matrix<T> LBFGS<T>::_g = Matrix<T>(256,256);
        template<class T> typename LBFGS<T>::RT LBFGS<T>::_rmse = 0;
        
        
    }}

#endif // _LBFGS_HPP_

