/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __ESTIMATE_SENSITIVITIES_RADIAL_VIBE_HPP__
#define __ESTIMATE_SENSITIVITIES_RADIAL_VIBE_HPP__

#include "ReconStrategy.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief Empty recon for test purposes
	 */
	class EstimateSensitivitiesRadialVIBE : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		EstimateSensitivitiesRadialVIBE  () : _cart_3rd_dim(false) {}
		
		/**
		 * @brief Default destructor
		 */
		virtual ~EstimateSensitivitiesRadialVIBE () {}
		
		
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

		Vector<size_t> _image_space_dims;
		bool _cart_3rd_dim;

	};

}
#endif /* __DUMMY_RECON_H__ */

