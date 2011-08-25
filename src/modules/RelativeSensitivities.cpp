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

	Matrix<cplx>*   data = m_cplx["data"];
	Matrix<cplx>*   shim = m_cplx["shim"];
	Matrix<cplx>*   rxm  = m_cplx["rxm"];
	Matrix<cplx>*   txm  = m_cplx["txm"];
	Matrix<double>* snro = m_real["snro"];
	Matrix<double>* b0   = m_real["b0"];

	printf ("  Processing map generation ...\n");

	// Squeeze matrix ---------------------------

	data->Squeeze();
	printf ("  Dimensions: %s \n", data->DimsToCString());
	// -----------------------------------------


	// Fourier transform -----------------------

	FTVolumes (data);
	// -----------------------------------------
	

	// Remove readout oversampling -------------

	RemoveOS (data);
	// -----------------------------------------


	// SVD calibration -------------------------
	
	shim  = new Matrix<cplx>   (data->Dim(4), 1);                                        // Shim coefficients
	txm   = new Matrix<cplx>   (data->Dim(0), data->Dim(1), data->Dim(2), data->Dim(4)); // TX maps
	rxm   = new Matrix<cplx>   (data->Dim(0), data->Dim(1), data->Dim(2), data->Dim(5)); // RX maps
	snro  = new Matrix<double> (data->Dim(0), data->Dim(1), data->Dim(2));               // SNR optimal image
	
	SVDCalibrate (data, rxm, txm, snro, shim, false);
   	// -----------------------------------------


	// B0 calculation --------------------------

	b0       = new Matrix<double> (data->Dim(0), data->Dim(1), data->Dim(2));
	float TE = 1.5e-3;

	B0Map (data, b0, TE);
	// -----------------------------------------

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
