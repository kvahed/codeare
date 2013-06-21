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

#include "matrix/io/IOContext.hpp"

template <class T> bool
cstest (Connector<T>* con) {

    using namespace codeare::matrix::io;

	Matrix<cxfl>  indata;
	Matrix<cxfl>  im_dc;
	Matrix<float> mask;
	Matrix<float> pdf;
	Matrix<cxfl>  pc;
	
    // Read data
    IOContext ic (con->GetElement("/config/data-in"), base, READ);
    indata = ic.Read<cxfl>(con->GetElement("/config/data-in/data"));
    mask   = ic.Read<float>(con->GetElement("/config/data-in/mask"));
    pdf    = ic.Read<float>(con->GetElement("/config/data-in/pdf"));
    pc     = ic.Read<cxfl>(con->GetElement("/config/data-in/pc"));
    ic.Close();

	if (con->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}
	
	// Outgoing -------------
	
	con->SetMatrix  ("data", indata); // Measurement data
	con->SetMatrix  ("pdf",  pdf);    // Sensitivities
	con->SetMatrix  ("mask", mask);   // Weights
	con->SetMatrix  ("pc",   pc);     // Phase correction
	// ---------------------
	
	con->Process (test);
	
	// Incoming -------------
	
	con->GetMatrix ("im_dc", im_dc);  // Recon output
	// ---------------------
	
	con->Finalise   (test);
	
    // Write images
    IOContext oc (con->GetElement("/config/data-out"), base, WRITE);
    oc.Write(im_dc, "image");
    oc.Close();
	
	return true;
	
}

