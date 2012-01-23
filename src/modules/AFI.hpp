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




	/**
	 * @brief              Noise prewhitening
	 * 
	 * @param  data    Data
	 * @param  ncov    Noise covariance matrix
	 * 
	 */
	void PreWhite (const Matrix<cxfl>& in, const Matrix<cxfl>& ncov, Matrix<cxfl>& od, const size_t& out = -1) {

		//assert (idx <= in.HDim());

		


	} 

		/*

imsz = size(data);

if imsz(channel_idx) == 1
    error([TAG 'Only single channel data?']);
end

% If just variances are given transform to diagonal matrix
if isvector(noise_covar)  
    noise_covar = diag(noise_covar);
end

if ~ndims(noise_covar) == 2 || (size(noise_covar,1) ~= size(noise_covar,2))
     error([TAG 'Invalid size of the covariance matrix! Has to be Nch x Nch!']);
end

if size(noise_covar,1) ~= imsz(channel_idx)
    error([TAG 'Invalid size of the covariance matrix! Has to be Nch x Nch!']);
end


%Make channel first index
old_order = 1:ndims(data);
new_order = [channel_idx old_order(old_order ~= channel_idx)]; 
data = permute(data, new_order); 
sz = size(data);

%Prewhite the data
R = (sqrtm(noise_covar));
whitened = R\data(:,:);
whitened = reshape(whitened,sz);

%Bring back to original shape
whitened = ipermute(whitened,new_order);
		*/


	void PhaseCombine (const Matrix<cxfl>& img, Matrix<float>& pc) {
		
		size_t nc = img.HDim();

		if (pc.IsZero())
			return;

	}



	void PhasePreset (const Matrix<cxfl>& afid, const bool& use_real, Matrix<float>& phase) {
		
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
							afid (x, y, z, 1) *= exp (cxfl (0.0, pc(x, y, z)));		

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

}
#endif /* __AFI_H__ */

