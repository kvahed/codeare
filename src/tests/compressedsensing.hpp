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

template <class T> bool
cstest (Connector<T>* con) {

	Matrix<cxfl>  indata;
	Matrix<cxfl>  im_dc;
	Matrix<float> mask;
	Matrix<float> pdf;
	Matrix<cxfl>  pc;
	
	std::string   cf  = std::string (base + std::string(config));
	std::string   df  = std::string (base + std::string(data));
	std::string   odf = std::string (base + std::string("/csout.mat"));

#ifdef HAVE_MAT_H	
	if (!(MXRead (indata, df, "data")))	return false;
	if (!(MXRead (pdf,    df, "pdf")))	pdf  = Matrix<float>(1);
	if (!(MXRead (mask,   df, "mask")))	mask = Matrix<float>(1);
	if (!(MXRead (pc,     df, "ph")))   pc   = Matrix<cxfl>(1);
#endif

	con->ReadConfig (cf.c_str());
	
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
	con->GetMatrix ("data",  indata); // Weighted FT of original input
	
	// ---------------------
	
	con->Finalise   (test);
	
#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen (odf.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}

	MXDump (indata, mf, "in");
	MXDump (im_dc,  mf, "out");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}

#endif
	
	return true;
	
}

