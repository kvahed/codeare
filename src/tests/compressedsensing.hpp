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

template <class T> bool
cstest (Connector<T>* rc) {

	Matrix<cxfl> indata;
	Matrix<cxfl> im_dc;
	Matrix<double> mask;
	Matrix<double> pdf;
	
	std::string cf  = std::string (base + std::string(config));
	std::string df  = std::string (base + std::string(data));
	std::string odf = std::string (base + std::string("/csout.mat"));

#ifdef HAVE_MAT_H	
	indata.MXRead (df, "data");
	pdf.MXRead (df, "pdf");
	mask.MXRead (df, "mask");
#endif

	rc->ReadConfig (cf.c_str());
	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}
	
	// Outgoing -------------
	
	rc->SetMatrix  ("data", indata); // Measurement data
	rc->SetMatrix  ("pdf",  pdf);  // Sensitivities
	rc->SetMatrix  ("mask", mask); // Weights
	
	// ---------------------
	
	rc->Process (test);
	
	// Incoming -------------
	
	rc->GetMatrix ("data", indata);  // Images
	rc->GetMatrix ("im_dc", im_dc);  // Images
	rc->GetMatrix ("orig", im_dc);  // Images
	
	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (odf.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}

	indata.MXDump (mf, "img");
	im_dc.MXDump (mf, "wvt");
	im_dc.MXDump (mf, "orig");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif
	
	return true;
	
}

