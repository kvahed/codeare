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
cgsensetest (RRClient::Connector<T>* rc) {
    
    using namespace codeare::matrix::io;

    // PARC or being scanned
    
    // Incoming
    Matrix<cxfl>  rawdata; // Measurement data O(Nkx,Nky,Nkz)
    Matrix<cxfl>  sens;    // Sensitivity maps O(Nx, Ny, Nz)
    Matrix<float> weights; // NUFFT weights
    Matrix<float> kspace;  // Kspace positions O(Nkx,Nky,Nkz)

    // Read data
    IOContext ic (rc->GetElement("/config/data-in"), base, READ);
    rawdata = ic.Read<cxfl>(rc->GetElement("/config/data-in/d"));
    kspace  = ic.Read<float>(rc->GetElement("/config/data-in/k"));
    weights = ic.Read<float>(rc->GetElement("/config/data-in/w"));
    sens	= ic.Read<cxfl>(rc->GetElement("/config/data-in/s"));
    ic.Close();

    rc->SetAttribute ("hans", 10.0);
    
    if (rc->Init (test) != codeare::OK) {
        printf ("Intialising failed ... bailing out!"); 
        return false;
    }
    
    // Send one time data
    rc->SetMatrix ("weights", weights); // Weights
    rc->SetMatrix ( "kspace", kspace);  // K-space
    rc->SetMatrix (   "sens", sens);    // Sensitivities
    rc->Prepare   (test);

    // Recon
    rc->SetMatrix (   "data", rawdata); // Measurement data
    rc->Process   (test);

    // Outgoing
    Matrix<cxfl>  cgimg;

    rc->GetMatrix (  "image", cgimg);  // Images
    rc->Finalise   (test);

    // Write images

	std::string ofname = base;
    ofname += "/";
	ofname += rc->GetElement("/config/data-out")->Attribute("fname");
	IOContext ofile = fopen (ofname, WRITE);
	fwrite (ofile, cgimg);
	fclose (ofile);

    return true;
    
}
