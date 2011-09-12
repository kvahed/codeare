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

	Matrix<cplx>*   data = m_cplx["meas"];
	Matrix<cplx>*   mask = m_cplx["mask"];

	printf ("  Processing map generation ...\n");
	ticks start = getticks();


	// Squeeze matrices ---------------------------

	data->Squeeze();
	mask->Squeeze();

	printf ("  Data:\n");
	printf ("    Dimensions: %s \n", data->DimsToCString());
	printf ("    Reolutions: %s \n", data->ResToCString());
	printf ("  Mask:\n");
	printf ("    Dimensions: %s \n", mask->DimsToCString());
	printf ("    Reolutions: %s \n", mask->ResToCString());
	// -----------------------------------------

	// Fourier transform -----------------------

	FTVolumes (data);
	RemoveOS (data);
	// -----------------------------------------

	// SVD calibration -------------------------

	Matrix<cplx>*   rxm;
	Matrix<cplx>*   txm;
	Matrix<cplx>*   shim;
	Matrix<double>* snro;

	AddCplx ("rxm",  rxm  = new Matrix<cplx>   (data->Dim(0), data->Dim(1), data->Dim(2), data->Dim(5)));
	AddCplx ("txm",  txm  = new Matrix<cplx>   (data->Dim(0), data->Dim(1), data->Dim(2), data->Dim(4)));
	AddCplx ("shim", shim = new Matrix<cplx>   (data->Dim(4), 1));
	AddReal ("snro", snro = new Matrix<double> (data->Dim(0), data->Dim(1), data->Dim(2)));

	SVDCalibrate (data, rxm, txm, snro, shim, false);

   	// -----------------------------------------

	// Do we have GRE for segmentation? --------

	FTVolumes (mask);
	RemoveOS  (mask);
	mask->SOS (mask->HDim());
	
	Matrix<double> bet(mask->Dim());

	Matrix<short>* bets;
	AddPixel ("bets", bets = new Matrix<short> (mask->Dim()));

	double tmp;
	for (int i = 0; i < mask->Size(); i++) {
		tmp = log(abs(mask->At(i)));
		bet.At(i) = (tmp < -2.5) ? 0.0 : tmp + 2.5;
	}
	
	SegmentBrain(&bet, bets);
	bets->Resample(0.5, LINEAR);

	// -----------------------------------------

	// B0 calculation --------------------------

	Matrix<double>* b0;
	AddReal ("b0", b0 = new Matrix<double> (data->Dim(0), data->Dim(1), data->Dim(2)));

	B0Map (data, b0, m_echo_shift);
	// -----------------------------------------

	for (int ch = 0; ch < txm->Dim(3); ch++)
		for (int i = 0; i < bets->Size(); i++)
			txm->At(ch*bets->Size() + i) *= (double)bets->At(i);
	
	for (int ch = 0; ch < rxm->Dim(3); ch++)
		for (int i = 0; i < bets->Size(); i++)
			rxm->At(ch*bets->Size() + i) *= (double)bets->At(i);
	
	for (int i = 0; i < bets->Size(); i++)
		b0->At(i) *= (double)bets->At(i);
	
	// Remove original data --------------------

    FreeCplx ("meas");
	// -----------------------------------------

	
	printf ("... done. Overall WTime: %.4f seconds.\n\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());
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
