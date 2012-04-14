/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "RelativeSensitivities.hpp"
#include "MITK.hpp"


using namespace RRStrategy;


RRSModule::error_code
RelativeSensitivities::Init        () {

    Attribute ("echo_shift", &m_echo_shift);
	printf ("  dTE[ms]: %.4f\n", m_echo_shift);
    Attribute ("cutoff",     &m_cutoff);
	printf ("  Cut off threshold: %.4f\n", m_cutoff);
    Attribute ("use_bet",    &m_use_bet);
	printf ("  Use FSL bet: %i\n", m_use_bet);
    Attribute ("log_mask",   &m_log_mask);
	printf ("  Log mask: %i\n", m_log_mask);
	Attribute ("weigh_maps", &m_weigh_maps);
	printf ("  Mask maps: %i\n", m_weigh_maps);
    
    return RRSModule::OK;

}


RRSModule::error_code
RelativeSensitivities::Process     () { 

    Matrix<cxfl>& data = GetCXFL("meas");
    Matrix<cxfl>& mask = GetCXFL("mask");

    printf ("  Processing map generation ...\n");
    ticks start = getticks();


    // Squeeze matrices ---------------------------

	data = squeeze(data);
    printf ("  Data:\n");
    printf ("    Dimensions: %s \n", DimsToCString(data));
    printf ("    Reolutions: %s \n", ResToCString(data));

	mask = squeeze(mask);
    printf ("  Mask:\n");
    printf ("    Dimensions: %s \n", DimsToCString(mask));
    printf ("    Reolutions: %s \n", ResToCString(mask));
    // -----------------------------------------

    // Fourier transform -----------------------

    FTVolumes (data);
    RemoveOS (data);
    // -----------------------------------------

    // SVD calibration -------------------------
    Matrix<cxfl>&   rxm  = AddMatrix ("rxm",  (Ptr<Matrix<cxfl> >)   NEW (Matrix<cxfl>   (data.Dim(0), data.Dim(1), data.Dim(2), data.Dim(5))));
    Matrix<cxfl>&   txm  = AddMatrix ("txm",  (Ptr<Matrix<cxfl> >)   NEW (Matrix<cxfl>   (data.Dim(0), data.Dim(1), data.Dim(2), data.Dim(4))));
    Matrix<cxfl>&   shim = AddMatrix ("shim", (Ptr<Matrix<cxfl> >)   NEW (Matrix<cxfl>   (data.Dim(4), 1)));
    Matrix<double>& snro = AddMatrix ("snro", (Ptr<Matrix<double> >) NEW (Matrix<double> (data.Dim(0), data.Dim(1), data.Dim(2))));

    SVDCalibrate (data, rxm, txm, snro, shim, false);

	// -----------------------------------------

    // Do we have GRE for segmentation? --------
	
    Matrix<short>& bets = AddMatrix ("bets", (Ptr<Matrix<short> >) NEW (Matrix<short> (mask.Dim())));
	
    if (m_use_bet == 1) { // Better test? / Replace with SNRO?
		
        FTVolumes (mask);
        RemoveOS  (mask);

		SOS (mask, ndims(mask));
        
        Matrix<double> bet(mask.Dim());
        double         tmp = 0.0;
		
        for (size_t i = 0; i < numel(mask); i++) {
            tmp = log(abs(mask[i]));
            bet[i] = (tmp < m_cutoff) ? 0.0 : tmp - m_cutoff;
        }
        
        SegmentBrain (bet, bets);
		bets = Resample (bets, 0.5, LINEAR);
		
    } else if (m_use_bet == 2) {
		
	    LogMask (snro, m_cutoff);
		bets = (Matrix<short>) snro;
		
    } else {
		
        for (size_t i = 0; i < bets.Size(); i++)
            bets[i] = 1;
		
    }
	
    // -----------------------------------------

    // B0 calculation --------------------------
	
    Matrix<double>& b0 = AddMatrix ("b0", (Ptr<Matrix<double> >) NEW (Matrix<double> (data.Dim(0), data.Dim(1), data.Dim(2))));
	
    B0Map (data, b0, m_echo_shift);
    // -----------------------------------------

    // Weighing with masks ---------------------

	if (m_weigh_maps) {

		for (int ch = 0; ch < size(txm,3); ch++)
			for (int i = 0; i < numel(bets); i++)
				txm[ch*numel(bets) + i] *= (double)bets[i];
		
		for (int ch = 0; ch < size(rxm,3); ch++)
			for (int i = 0; i < numel(bets); i++)
				rxm[ch*numel(bets) + i] *= (double)bets[i];
		
		for (int i = 0; i < numel(bets); i++)
			b0[i] *= (double)bets[i];

	}
	// -----------------------------------------

    // Remove original data from RAM -----------

    FreeCXFL ("meas");
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
