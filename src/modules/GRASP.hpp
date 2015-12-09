/*
 *  codeare Copyright (C) 2007-2012 Kaveh Vahedipour
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

#ifndef __GRASP_HPP__
#define __GRASP_HPP__

#include "ReconStrategy.hpp"
#include "Algos.hpp"
#include "DFT.hpp"
#include "CS_XSENSE.hpp"
#include "DWT.hpp"
#include "TVOP.hpp"
#include "CX.hpp"
#include "linalg/Lapack.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {


    /**
     * @brief CS reconstruction based on Sparse MRI v0.2 by Michael Lustig
     */
    class GRASP : public ReconStrategy {
            
        
    public:
        
        /**
         * @brief Default constructor
         */
        GRASP  ();
        
        
        /**
         * @brief Default destructor
         */
        virtual ~GRASP ();
        
        
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
        
        
    private:
        
        Params ft_params;
        size_t _tf;
        float _ta;

    };
    
}
#endif /* __GRASP_H__ */

