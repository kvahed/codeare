/*
 *  codeare Copyright (C) 2007-2012 Kaveh Vahedipour
 *                                  Daniel Joergens
 *                                  Forschungszentrum Juelich, Germany
 *
 *  Stolen ;) from sparse MRI 0.2 by Michael Lustig
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#ifndef __COMPRESSED_SENSING_HPP__
#define __COMPRESSED_SENSING_HPP__

#include "ReconStrategy.hpp"
#include "Algos.hpp"
#include "FT.hpp"
#include "DWT.hpp"
#include "TVOP.hpp"
#include "CX.hpp"
#include "linalg/Lapack.hpp"

#ifdef _MSC_VER
std::string ofstr = "    %02Iu - nrms: %1.7f, l-search: %d, ";
#else
std::string ofstr = "    %02zu - nrms: %1.7f, l-search: %d, ";
#endif
//#include <pthread.h>
/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {


    /**
     * @brief CG parameters
     */
    struct CSParam {
        
        int    lsiter;
        int    lsto ;
        int    cgiter;
        
        double pnorm;
        double tvw;
        double xfmw;
        double cgconv;
        double l1;
        double lsa;
        double lsb;

        DWT<cxfl>* dwt;
        FT<float>*  ft;
        TVOP*      tvt;
        
    };

    
    /**
     * @brief CS reconstruction based on Sparse MRI v0.2 by Michael Lustig
     */
    class CompressedSensing : public ReconStrategy {
            
        
    public:
        
        /**
         * @brief Default constructor
         */
        CompressedSensing  ();
        
        
        /**
         * @brief Default destructor
         */
        virtual 
        ~CompressedSensing ();
        
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code
        Process ();
        
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code 
        Init ();
        
        /**
         * @brief Do nothing
         */
        virtual codeare::error_code
        Prepare ();
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code
        Finalise ();
        
        
        
    private:
        
        int            m_dim;    /**< Image recon dim */
        int            m_N[3];   /**< Data side lengths */
        int            m_csiter; /**< # global iterations */
        Vector<size_t> m_image_size;
        int            m_test_case;

        double m_noise;

        CSParam        m_csparam;
        int            m_wf;
        int            m_wm;
        int            m_verbose;
        
        int            m_ft_type;

    };
    
    
    inline static float 
    Obj (const Matrix<cxfl>& ffdbx, const Matrix<cxfl>& ffdbg, 
         const Matrix<cxfl>&  data, const float             t) {
        
        Matrix<cxfl> om = ffdbx;

        if (t > 0.0) 
            om += t * ffdbg;
        om -= data;

        return real(om.dotc(om));

    }
    
    
    inline static float 
    ObjTV (const Matrix<cxfl>& ttdbx, const Matrix<cxfl>& ttdbg, 
           const float             t, const CSParam&        cgp) {
        
        float o = 0.0, p = 0.5*cgp.pnorm;
        Matrix<cxfl> om = ttdbx;
        
        if (t > 0.0)
            om += t * ttdbg;
        om *= conj(om);
        om += cgp.l1;
        om ^= p;
        
        for (size_t i = 0; i < om.Size(); i++) 
            o += real(om[i]);
        
        return cgp.tvw * o;
        
    }
    
    
    /**
     *
     */
    inline static float 
    ObjXFM (const Matrix<cxfl>& x, const Matrix<cxfl>& g, 
            const float         t, const CSParam&    cgp) {
        
        float o = 0.0, p = 0.5*cgp.pnorm;
        Matrix<cxfl> om = x;

        if (t > 0.0)
            om += t * g;
        om *= conj(om);
        om += cgp.l1;
        om ^= p;
        
        for (size_t i = 0; i < om.Size(); i++) 
            o += om[i].real();
        
        return cgp.xfmw * o;
        
    } 
    
    
    static float 
    Objective (const Matrix<cxfl>& ffdbx, const Matrix<cxfl>& ffdbg, 
               const Matrix<cxfl>& ttdbx, const Matrix<cxfl>& ttdbg, 
               const Matrix<cxfl>&     x, const Matrix<cxfl>&     g, 
               const Matrix<cxfl>&  data, const float             t, 
                     float&         rmse, const CSParam&        cgp) {
        
        float obj = Obj (ffdbx, ffdbg, data, t);
        
        rmse = sqrt(obj/(float)nnz(data));
        if (cgp.tvw)
            obj += ObjTV (ttdbx, ttdbg, t, cgp);
        if (cgp.xfmw)
            obj += ObjXFM (x, g, t, cgp);
        
        return obj;

    }
    
    
    /**
     * @brief Compute gradient of the data consistency
     */
    static Matrix<cxfl> 
    GradObj (const Matrix<cxfl>& x, const Matrix<cxfl>& wx, 
             const Matrix<cxfl>& data, const CSParam& cgp) {
        
        FT<float>& ft = *(cgp.ft);
        DWT<cxfl>& dwt = *(cgp.dwt);

        return (2.0 * (dwt * (ft ->* ((ft * wx) - data))));
        
    }

    
    /**
     * @brief Compute gradient of L1-transform operator
     *
     * @param  x   X
     * @param  cgp CG parameters
     * @return     The gradient
     */
    template <class T> inline static Matrix<T> 
    GradXFM   (const Matrix<T>& x, const CSParam& cgp) {
        
        float pn = 0.5*cgp.pnorm-1.0, l1 = cgp.l1, xfm = cgp.xfmw;
        return xfm * (x * ((x * conj(x) + l1) ^ pn));
        
    }
    
    
    /**
     * @brief Compute gradient of the total variation operator
     *
     * @param  x   Image space original
     * @param  wx  Image space perturbance
     * @param  cgp Parameters
     */
    Matrix<cxfl> 
    GradTV    (const Matrix<cxfl>& x, const Matrix<cxfl>& wx, const CSParam& cgp) {
        
        DWT<cxfl>&  dwt = *cgp.dwt;
        TVOP& tvt = *cgp.tvt;
        float p   = 0.5*cgp.pnorm-1.0;
        Matrix<cxfl> dx, g;
        
        dx = tvt * wx;
        g  = dx * conj(dx);
        g += cgp.l1;
        g ^= p;
        g *= dx;
        g *= cgp.pnorm;
        g  = dwt * (tvt->*g);
        
        return (cgp.tvw * g);
        
    }
    
    
    Matrix<cxfl> 
    Gradient (const Matrix<cxfl>& x, const Matrix<cxfl>& wx, const Matrix<cxfl>& data, const CSParam& cgp) {

        Matrix<cxfl> g = GradObj (x, wx, data, cgp);
        return g += (cgp.xfmw) ? GradXFM (x, cgp) : GradTV  (x, wx, cgp);

    } 
    

    void NLCG (Matrix<cxfl>& x, const Matrix<cxfl>& data, const CSParam& cgp) {

        
        float     t0  = 1.0, t = 1.0, z = 0., xn = norm(x), rmse, bk, f0, f1, dxn;
        
        Matrix<cxfl> g0, g1, dx, ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx;
        
        DWT<cxfl>& dwt = *cgp.dwt;
        FT<float>& ft  = *cgp.ft;
        TVOP&      tvt = *cgp.tvt;
        
        wx  = dwt->*x;
        g0 = Gradient (x, wx, data, cgp);
        dx = -g0;
        wdx = dwt->*dx;

        for (size_t k = 0; k < (size_t)cgp.cgiter; k++) {
            
            t = t0;
            
            ffdbx = ft * wx;
            ffdbg = ft * wdx;
            
            if (cgp.tvw) {
                ttdbx = tvt * wx;
                ttdbg = tvt * wdx;
            }
            
            f0 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, z, rmse, cgp);
            
            int i = 0;
            while (i < cgp.lsiter) {
                
                t *= cgp.lsb;
                f1 = Objective(ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, t, rmse, cgp);
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
            x  += (dx * t);
            wx  = dwt->*x;
                    
            // CG computation 
            g1  =  Gradient (x, wx, data, cgp);
            bk  =  real(g1.dotc(g1)) / real(g0.dotc(g0));
            g0  =  g1;
            dx  = -g1 + dx * bk;
            wdx =  dwt->*dx;
            dxn =  norm(dx)/xn;
            
            printf ("dxnrm: %0.4f\n", dxn);
            if (dxn < cgp.cgconv) 
                break;
            
        } 

    }
    
}
#endif /* __COMPRESSED_SENSING_H__ */

