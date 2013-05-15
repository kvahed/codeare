/*
 *  codeare Copyright (C) 2010-2012 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
 *  
 *  Original code stolen from Martijn Cloos
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

#ifndef __KT_POINTS_HPP__
#define __KT_POINTS_HPP__

#include "ReconStrategy.hpp"
#include "Toolbox.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IO.hpp"
#include "linalg/Lapack.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

    /**
     * @brief Spatial domain method for RF pulse generation with variable exchange method<br/>
     *        Cloos MA, Boulant N, Luong M, Ferrand G, Giacomini E, Le Bihan D, Amadon A.<br/>
     *        kT-points: short three-dimensional tailored RF pulses for flip-angle homogenization over an extended volume. MRM:2012; 67(1), 72-80.
     */
    class KTPoints : public ReconStrategy {
        
        
    public:
        
        /**
         * @brief Default constructor
         */
        KTPoints  ();
        
        /**
         * @brief Default destructor
         */
        virtual 
        ~KTPoints ();
        
        
        /**
         * @brief Process
         */
        virtual error_code
        Process ();

        
        /**
         * @brief Initialise
         */
        virtual error_code
        Init ();
        

        /**
         * @brief Finalise
         */
        virtual error_code
        Finalise ();



    private:
 
        int           nc;       /**< @brief Transmit channels       */
        Matrix<short> m_pd;       /**< @brief Pulse durations         */
        int           m_gd;       /**< @brief Gradient durations      */
        int           ns;       /**< @brief # Spatial positions     */
        int           nk;       /**< @brief # kt-points             */
        int           m_maxiter;  /**< @brief # Variable exchange method iterations */
        int           m_verbose;  /**< @brief Verbose output. All intermediate results. */
        int           m_breakearly;  /**< @brief Break search with first diverging step */

        double        m_lambda;   /**< @brief Tikhonov parameter      */
        double        m_rflim;    /**< @brief Maximum rf amplitude    */
        double        m_conv;     /**< @brief Convergence criterium   */
        
        float*        m_max_rf;   /**< @brief Maximum reached RF amps */

        std::string   m_orient;   /**< @brief Orientation             */ 
        std::string   m_ptxfname; /**< @brief PTX file name           */

    };

}
#endif /* __KT_POINTS_HPP__ */


