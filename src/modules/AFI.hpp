/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Daniel Brenner
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

#ifndef __AFI_HPP__
#define __AFI_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief AFI reconstruction
	 */
	class AFI : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		AFI  ();
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~AFI ();
		
		
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

		bool   m_use_real;
		bool   m_retain_phase;
		
	};



	void PhaseCombine (const Matrix<cplx>& img, Matrix<float>& pc) {
		
		size_t nc = img.HDim();

		if (pc.IsZero())
			return;

	}


	void PhasePreset (const Matrix<cplx>& afid, const bool& use_real, Matrix<float>& phase) {
		
		if (afid.Dim(4) > 1) {
			
			if (use_real) {
				
				for (size_t z = 0; z < afid.Dim(2); z++)
					for (size_t y = 0; y < afid.Dim(1); y++)
						for (size_t x = 0; x < afid.Dim(0); x++)
							phase (x, y, z) = arg(afid(x, y, z, 0, 0));

				Matrix<float> pc (afid.Dim(0), afid.Dim(1), afid.Dim(2));
				PhaseCombine (afid, pc);

				afid.SOS(4);
				
				for (size_t z = 0; z < afid.Dim(2); z++)
					for (size_t y = 0; y < afid.Dim(1); y++)
						for (size_t x = 0; x < afid.Dim(0); x++)
							afid (x, y, z, 1) *= exp (cplx (0.0, pc(x, y, z)));		

			} else {
				
				afid.SOS(4);
				
			}

		} else {

			for (size_t z = 0; z < afid.Dim(2); z++)
				for (size_t y = 0; y < afid.Dim(1); y++)
					for (size_t x = 0; x < afid.Dim(0); x++)
						phase (x,y,z) = arg(afid(x,y,z,0,0));

		}
		
	}


			/*
    % Multi coil data set
    if ndims(images) > 4   % Means we have a multicoil data set
       if use_real
             phase               = angle(images(:,:,:,1,1));
           ph              = phase_combine(images);
           images          = ssq_rec(images);
           images(:,:,:,2) = images(:,:,:,2).*exp(1i*ph);
       else
           images          = ssq_rec(images);
       end
       
     
    else
       phase               = angle(images(:,:,:,1)); 
    end
	*/


}
#endif /* __AFI_H__ */

