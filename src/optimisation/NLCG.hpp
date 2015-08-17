/*
 * NLCG.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _NLCG_HPP_
#define _NLCG_HPP_

#include <NonLinear.hpp>

namespace codeare {
    namespace optimisation {
        
template<class T> class NLCG : public NonLinear<T> {

public:
	NLCG (const size_t& iterations) : NonLinear<T>::NonLinear (iterations)  {};
	virtual ~NLCG() {};
    
	inline Matrix<T> Minimize(/*const Operator<T>& A, */const MatrixType<T>& x, const MatrixType<T>& data) {
/*        
        typename TypeTraits<T>::RT t0  = 1.0, t = 1.0, z = 0., xn = norm(x), rmse, bk, f0, f1, dxn;
        Matrix<T> dx, g0, g1;
        
        g0  = df (x, data, cgp);
        dx  = -g0;
        
        for (size_t k = 0; k < (size_t)cgp.cgiter; k++) {

            Update(dx, cgp);
            
            t = t0;
            
            f0 = f (x, dx, data, z, rmse, cgp);
            
            int i = 0;
            while (i < cgp.lsiter) {
                t *= cgp.lsb;
                f1 = f (x, dx, data, t, rmse, cgp);
                if (f1 <= f0 - (cgp.lsa * t * abs(g0.dotc(dx))))
                    break;
                ++i;
            }
            
            printf (ofstr.c_str(), k, rmse, i); fflush (stdout);
            
            if (i == cgp.lsiter) {
                printf ("Reached max line search, exiting... \n");
                return;
            }
            
            if      (i > 2) t0 *= cgp.lsb;
            else if (i < 1) t0 /= cgp.lsb;
            
            // Update image
            x  += dx * t;
            
            // CG computation
            g1  =  df (x, data, cgp);
            bk  =  real(g1.dotc(g1)) / real(g0.dotc(g0));
            g0  =  g1;
            dx  = -g1 + dx * bk;
            dxn =  norm(dx)/xn;
            
            printf ("dxnrm: %0.4f\n", dxn);
            if (dxn < cgp.cgconv)
                break;
            
        }
*/        
    }
        
private:
    
/*    typename TypeTraits<T>::RT t0  = 1.0, t = 1.0, z = 0., xn, rmse, bk, f0, f1, dxn;
    Matrix<cxfl> g0, g1, dx, ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx;
    const DWT<cxfl>&  dwt = *cgp.dwt;
    const FT<cxfl>&   ft  = *cgp.ft;
    const TVOP<cxfl>& tvt = *cgp.tvt;
*/  

};
        
    }}

#endif // _NLCG_HPP_
