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
RelativeSensitivities::Init        () {

	Attribute ("echo_shift", &m_echo_shift);

	return RRSModule::OK;

}


RRSModule::error_code
RelativeSensitivities::Process     () { 

	// Introductive output --------------------

	printf ("  Processing map generation ...\n");
	m_raw.Squeeze();

	printf ("  Dimensions: ");
	for (int i = 0; i < INVALID_DIM; i++)
		if (m_raw.Dim(i) > 1)
			printf (" %i",  m_raw.Dim(i));
	printf ("\n");
	// ----------------------------------------


	// Fourier transform ----------------------

	FTVolumes (&m_raw);
	// -----------------------------------------
	

	// Remove readout oversampling -------------

	RemoveOS (&m_raw);
	// -----------------------------------------


	// SVD calibration -------------------------
	
	Matrix<raw>    shim (m_raw.Dim(4), 1);                                        // Shim coefficients
	Matrix<raw>    txm  (m_raw.Dim(0), m_raw.Dim(1), m_raw.Dim(2), m_raw.Dim(4)); // TX maps
	Matrix<raw>    rxm  (m_raw.Dim(0), m_raw.Dim(1), m_raw.Dim(2), m_raw.Dim(5)); // RX maps
	Matrix<double> snro (m_raw.Dim(0), m_raw.Dim(1), m_raw.Dim(2));               // SNR optimal image
	
	SVDCalibrate (&m_raw, &rxm, &txm, &snro, &shim, false);
   	// -----------------------------------------


	// B0 calculation --------------------------

	Matrix<double> b0 (m_raw.Dim(0), m_raw.Dim(1), m_raw.Dim(2));
	float          TE = 1.5e-3;

	B0Map (&m_raw, &b0, TE);
	// -----------------------------------------

	m_raw     = txm;
	m_rhelper = rxm;
	m_helper  = snro;
	m_kspace  = b0;

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
