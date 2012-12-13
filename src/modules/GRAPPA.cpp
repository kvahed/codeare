/*
 *  codeare Copyright (C) 2010-2011 Kaveh Vahedipour
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

#include "GRAPPA.hpp"
//#include "CGRAPPA.hpp"
using namespace RRStrategy;

RRSModule::error_code
GRAPPA::Init () {
	
	printf ("Intialising %s ...\n", Name());
	
	// Intialising ---------------------
	
	m_nc       = 0;
	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;
	m_acsinc   = 0;
	
	for (int i = 0; i < 3; i++) {
		m_R[i] = 1; m_acs_dim[i] = 1; m_data_dim[i] = 1;
	}

	// ---------------------------------
	
	// Dimensions  ---------------------
	
	Attribute("dim",       &m_dim);
    printf ("  # dimensions: %i \n", m_dim);
	
	if (m_dim < 2 || m_dim > 3) {
		printf ("%s only supports 2-3 dimensions. %i was specified.\n", Name(), m_dim);
		return UNSUPPORTED_DIMENSION;
	}
	// ---------------------------------
	
	// Coils and accelerqtion ----------
	
	Attribute("nc",       &m_nc);
    printf ("  # channels: %i \n", m_nc);

	if (m_nc == 0) {
		printf ("Initialising %s with %i channels? Check configuration! FAILED: Bailing out!\n", Name(), m_nc);
		return CGSENSE_ZERO_CHANNELS;
	}
	
	int tmp = 0;

	// Acceleration factors
	Attribute ("rx",        &tmp); m_R[0]        = tmp;
	Attribute ("ry",        &tmp); m_R[1]        = tmp;
	Attribute ("rz",        &tmp); m_R[2]        = tmp;
    printf ("  acceleration factors: %.1f,%.1f,%.1f \n", m_R[0], m_R[1], m_R[2]);

	// Kernel dimensions
	Attribute ("kernx",     &tmp); m_kern_dim[0] = tmp;
	Attribute ("kerny",     &tmp); m_kern_dim[1] = tmp;
	Attribute ("kernz",     &tmp); m_kern_dim[2] = tmp;
    printf ("  kernel dimensions: %i,%i,%i \n", m_kern_dim[0], m_kern_dim[1], m_kern_dim[2]);

	// ACS dimensions
	Attribute ("acs_x",     &tmp); m_acs_dim[0]  = tmp;
	Attribute ("acs_y",     &tmp); m_acs_dim[1]  = tmp;
	Attribute ("acs_z",     &tmp); m_acs_dim[2]  = tmp;
    printf ("  acs dimesions: %i,%i,%i \n", m_acs_dim[0], m_acs_dim[1], m_acs_dim[2]);

	// Incoming data dimensions 
	// (Outgoing: m_data_dim[i] * m_R[i])
	Attribute ("nx",        &tmp); m_data_dim[0] = tmp;
	Attribute ("ny",        &tmp); m_data_dim[1] = tmp;
	Attribute ("nz",        &tmp); m_data_dim[2] = tmp;
    printf ("  measured dims: %i,%i,%.i \n", m_data_dim[0], m_data_dim[1], m_data_dim[2]);

	m_d[0] = floor ((float) m_kern_dim[0] / 2);
	m_d[1] =       ((float) m_kern_dim[1] / 2 - 1) * m_R[1]; 
    printf ("  steps: %.i,%.i \n", m_d[0], m_d[1]);

	// ---------------------------------

	// Test and verbosity --------------

	Attribute ("testcase", &m_testcase);
	Attribute ("verbose",  &m_verbose); 
	Attribute ("noise",    &m_noise);
	// ---------------------------------

	m_initialised = true;

	printf ("... done.\n\n");

	return RRSModule::OK;

}


RRSModule::error_code
GRAPPA::Prepare     () { 

	printf ("  Preparing %s ...\n", Name());

	Matrix<cxfl>& acs  = Get<cxfl>("acs");
	/*ComputeWeights (m_nc, m_acs_dim, m_kern_dim, m_d, m_R, acs, m_weights);*/

	
	
	printf ("... done.\n\n");
	return RRSModule::OK;
	

}

// On entry ------------------------
//
// m_raw:     Measured k-spaces      O (Nc x NKx/RX x NKy/RY x NKz/RZ)
// m_rhelper: ACS scans if external  O (Nc x NACSX  x NACSY  x NACSZ)
// ---------------------------------

// On exit -------------------------
//
// m_raw:     Reconstructed k-spaces O (Nc x NKx    x NKy    x NKz)
// ---------------------------------

RRSModule::error_code
GRAPPA::Process     () { 

	ticks cgstart = getticks();

	printf ("  Processing GRAPPA ...\n");

	Matrix<cxfl>& acs  = Get<cxfl>("acs");
	Matrix<cxfl>& data = Get<cxfl>("data");

	printf ("  data dims: %s\n", DimsToCString(data));

	printf ("... done.. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	return RRSModule::OK;

}



RRSModule::error_code
GRAPPA::Finalise () {

	return RRSModule::OK;

}



// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new GRAPPA;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}
