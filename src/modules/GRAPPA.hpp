/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#ifndef __GRAPPA_HPP__
#define __GRAPPA_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief GRAPPA PPI reconstruction
	 */
	class GRAPPA : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		GRAPPA  () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~GRAPPA () {};
		
		
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
		Prepare ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Finalise ();


	private: 
		
		
		int           m_dim;      /**< Reco dimension {2,3}                   */
		int           m_verbose;  /**< Verbose output                         */
		int           m_nc;       /**< # of coils                             */
		int           m_testcase; /**< Test case                              */
		int           m_noise;    /**< Add noise to signal (for testing)      */
		int           m_acsinc;   /**< ACS scans included                     */
		int           m_ns;       /**< # of source points                     */
		int           m_nt;       /**< # of target points                     */
		int           m_nr;       /**< # of kernel reps in ACS                */
		int           m_data_dim[3]; /**< Kernel dimensions                   */
		int           m_acs_dim[3];  /**< Kernel dimensions                   */
		int           m_kern_dim[3]; /**< Kernel dimensions                   */
		int           m_d[3];        /**< Step size                           */

		float         m_R[3];        /**< Acceleration factor in each dim     */

		Matrix<cplx>  m_weights;  /**< Weights                                */
		



		void ComputeWeights (const int* acs_dim, const int* kern_dim, const int* d, const float* R, const Matrix<cplx>* acs, Matrix<cplx>& weights) {

			ticks       tic     = getticks();
			printf ("  (Re-)Computing weights \n");
			
			// # Kernel repititios in ACS
			int nr = (acs_dim[1] - (kern_dim[1]-1) * R[1]) * (m_acs_dim[0] - (m_kern_dim[0] - 1));
			int ns = m_nc * kern_dim[0] * kern_dim[1]; // # Source points
			int nt = m_nc * (m_R[1]-1);      // # Target points
			
			Matrix<cplx> s (ns, nr);          // Source 
			Matrix<cplx> t (nt, nr);          // Target

			int c = 0;

			//yind=1:nyacs-(srcy-1)*af
			/*
			for (int i = d[0]; i < acs_dim[0] - d[0]; i++)
				for (int j = 0; j < acs_dim[1] - (kern_dim[1]-1) * R[1]; j++, c++) 
					for (int c = 0; c < acs_dim[2], c++) {
					    s(c+j*acs_dim[2]+i*acs_dim[1] - (kern_dim[1]-1) * R[1]) = acs->At(i,j,c)
							//memcpy (&s[c*ns+i*], &acs)
							//memcpy (&s[0], &acs[0], s.Size() * sizeof(cplx));
							//src(:,cnt) = reshape(acs(:,yind:af:yind+(srcy-1)*af,xind-dx:xind+dx),nc*srcy*srcx,1);                   
							//memcpy (&t[0], &acs[0], t.Size() * sizeof(cplx));
							//trg(:,cnt) = reshape(acs(:,yind+1+dy:yind+dy+af-1,xind),nc*(af-1),1);
							}*/
			/*
			  s = s.Pinv();
			  weights = t->*s;
			*/
			printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
			
		}
		
	};
	
}
#endif /* __GRAPPA_H__ */

