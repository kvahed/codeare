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
#include "DFT.hpp"
#include "DWT.hpp"
#include "TVOP.hpp"
#include "CX.hpp"
#include "linalg/Lapack.hpp"

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
        FT<cxfl>*  ft;
        TVOP<cxfl>* tvt;
        
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
        virtual ~CompressedSensing ();
        
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code Process ();
        
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code Init ();
        
        /**
         * @brief Do nothing
         */
        virtual codeare::error_code Prepare ();
        
        /**
         * @brief Do nothing 
         */
        virtual codeare::error_code Finalise ();
        
        template<class T> typename TypeTraits<T>::RT obj (const Matrix<T>& x, const Matrix<T>& dx,
                            const typename TypeTraits<T>::RT& t, typename TypeTraits<T>::RT& rmse);
        template<class T> Matrix<T> df (const Matrix<T>& x);
        template<class T> void Update (const Matrix<T>& dx);
        template<class T> void NLCG (Matrix<T>& x);
        
        
    private:
        
        int            m_dim;    /**< Image recon dim */
        int            m_N[3];   /**< Data side lengths */
        int            m_csiter; /**< # global iterations */
        Vector<size_t> m_image_size;
        int            m_test_case;
        float          m_ndnz;

        double m_noise;

        CSParam        m_csparam;
        int            m_wf;
        int            m_wm;
        int            m_verbose;
        
        int            m_ft_type;

        Matrix<cxfl> data;

    };
    
}
#endif /* __COMPRESSED_SENSING_H__ */

