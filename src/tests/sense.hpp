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
sensetest (RRClient::Connector<T>* rc) {
	
	// Incoming
	Matrix<cxfl> fimgs;   // Folded images O (Nc, RO, PE/AF, PE2)
	Matrix<cxfl> smaps;   // Sensitivity maps O (Nc, RO, PE, PE2)
	
	// Outgoing
	Matrix<cxfl> ufimg;   // Reconstructed unaliased O (RO, PE, PE2)

	// Compute g-factor maps?
	bool         compgfm = false;
	rc->Attribute("compgfm", &compgfm);

	std::string    odf = std::string (base + std::string("/images.mat")); // Binary Ouput (images etc)
	
	if (!Read (fimgs, rc->GetElement("/config/data/in"), base) ||
		!Read (smaps, rc->GetElement("/config/data/sm"), base))
		return false;
	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}
	
	// Prepare -------------
	
	rc->SetMatrix ("smaps", smaps); // Sensitivities
	rc->SetMatrix ("fimgs", fimgs); // Weights
	
	rc->Prepare   (test);
	
	// Process -------------
	
	rc->Process   (test);
	
	// Receive -------------
	
	rc->GetMatrix ("image", ufimg);  // Images
	
	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (odf.c_str(), "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}
	
	MXDump     (ufimg, mf, "ufimg");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif

	return true;
	
}
