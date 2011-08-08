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

#include "RelativeSensitivities.hpp"

using namespace RRStrategy;

RRSModule::error_code
RelativeSensitivities::Process     () { 

	std::vector<long> rdim;
	short             found = 0;

	// Introductive output --------------------
	printf ("  Processing map generation ...\n");
	printf ("  Dimensions: ");
	for (int i = 0; i < INVALID_DIM; i++)
		if (m_raw.Dim(i) > 1) {
			rdim.push_back(m_raw.Dim(i));
			printf (" %i",  m_raw.Dim(i));
			found++;
		}
	printf ("\n");
	// ----------------------------------------

	ticks start = getticks();
	long  imsize  = rdim.at(0) * rdim.at(1) * rdim.at(2);
	int   vols    = m_raw.Size() / (imsize);

	// Fourier transform ----------------------
	printf ("  Fourier transforming %i volumes of %li ...\n", vols, imsize);
	
	int threads;

#pragma omp parallel default (shared) 

	{
		threads = omp_get_num_threads();
	}
	
	Matrix<raw>* mr = new Matrix<raw>[threads];

	for (int i = 0; i < threads; i++) {
		mr[i] = Matrix<raw>(rdim.at(0),rdim.at(1),rdim.at(2));
		//		mr[i].Dim(0) = rdim.at(0);
		//mr[i].Dim(1) = rdim.at(1);
		//mr[i].Dim(2) = rdim.at(2);
		//mr[i].Reset();
	}	

#pragma omp parallel default (shared)
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = vols / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < vols; i++) {
			memcpy (&mr[tid][0], &m_raw[i*imsize], imsize * sizeof(raw));
			mr[tid] = mr[tid].fftshift();
			mr[tid] = mr[tid].ifft();
			mr[tid] = mr[tid].ifftshift();
			memcpy (&m_raw[i*imsize], &mr[tid][0], imsize * sizeof(raw));
		}

	}
	
	for (int i = 0; i < threads; i++)
		mr[i].~Matrix<raw>();

	//delete mr;
	printf ("  ... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / ClockRate());
	// -----------------------------------------

	return RRSModule::OK;

}


RRSModule::error_code
RelativeSensitivities::Init        () {

	Attribute ("echo_shift", &m_echo_shift);

	return RRSModule::OK;

}


RRSModule::error_code
RelativeSensitivities::Finalise    () {

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new RelativeSensitivities;

}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

	delete p;

}
