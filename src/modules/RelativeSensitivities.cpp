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

	// Introductive output --------------------
	printf ("  Processing map generation ...\n");
	m_raw.Squeeze();

	printf ("  Dimensions: ");
	for (int i = 0; i < INVALID_DIM; i++)
		printf (" %i",  m_raw.Dim(i));
	printf ("\n");
	
	// ----------------------------------------

	ticks start  = getticks();
	long  imsize = m_raw.Dim(0) * m_raw.Dim(1) * m_raw.Dim(2);
	int   vols   = m_raw.Size() / (imsize);

	// Fourier transform ----------------------
	printf ("  Fourier transforming %i volumes of %lix%lix%li ...\n", vols, m_raw.Dim(0), m_raw.Dim(1), m_raw.Dim(2));
	
	int threads  = 0;

#pragma omp parallel default (shared) 
	{
		threads  = omp_get_num_threads();
	}
	
	/*
	Matrix<raw> mr[threads];
	
	for (int i = 0; i < threads; i++)
		mr[i] = Matrix<raw>(m_raw.Dim(0),m_raw.Dim(1),m_raw.Dim(2));
	*/

	Matrix<raw>* mr;

#pragma omp parallel default (shared) private (mr)
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = vols / omp_get_num_threads();

		mr = new Matrix<raw>(m_raw.Dim(0),m_raw.Dim(1),m_raw.Dim(2));
		
#pragma omp for schedule (dynamic, chunk)
		
		/*		for (int i = 0; i < vols; i++) {
			memcpy (&mr[tid][0], &m_raw[i*imsize], imsize * sizeof(raw));
			mr[tid] = mr[tid].FFTShift();
			mr[tid] = mr[tid].IFFT();
			mr[tid] = mr[tid].IFFTShift();
			memcpy (&m_raw[i*imsize], &mr[tid][0], imsize * sizeof(raw));
			}*/

		for (int i = 0; i < vols; i++) {
			memcpy (&mr->At(0), &m_raw[i*imsize], imsize * sizeof(raw));
			mr[tid] = mr->FFTShift();
			mr[tid] = mr->IFFT();
			mr[tid] = mr->IFFTShift();
			memcpy (&m_raw[i*imsize], &mr->At(0), imsize * sizeof(raw));
		}
		
		delete mr;
		
	}
	
	/*for (int i = 0; i < threads; i++)
	  mr[i].~Matrix<raw>();*/
	
	m_raw = m_raw.SOS(0);

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
