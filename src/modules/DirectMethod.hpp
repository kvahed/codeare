/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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

#ifndef __DIRECT_METHOD_HPP__
#define __DIRECT_METHOD_HPP__

#include "ReconStrategy.hpp"
#include "SimulationContext.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

    /**
     * @brief Empty recon for test purposes
     */
    class DirectMethod : public ReconStrategy {
        
        
    public:
        

        /**
         * @brief Default constructor
         */
        DirectMethod  () :
            m_np (1),
            m_verbose (true),
            m_dt (1.0),
            m_mode (0)
        {};
        

        /**
         * @brief Default destructor
         */
        virtual 
        ~DirectMethod () {

            this->Finalise();

        }
        
        
        /**
         * @brief Do nothing 
         */
        virtual RRSModule::error_code
        Process ();
        

        /**
         * @brief Do nothing 
         */
        virtual RRSModule::error_code
        Init ();
        
        /**
         * @brief Do nothing 
         */
        virtual RRSModule::error_code
        Finalise ();
        

    private: 
        
        double             m_dt;       /*!< @brief Simulation time steps                                        */
        bool               m_verbose;  /*!< @brief Verbose (Store magnetisation for every dt. HANDLE WITH CARE) */  // ignored
        int                m_np;       /*!< @brief Number of OMP threads */
        int                m_mode;     /*!< Single run: 0, Iterative: 1 */
		double             m_cgeps;
		int                m_cgiter;
		bool               m_cb0;  
		float              m_lambda;

    };


}
#endif /* __DIRECT_METHOD_H__ */

