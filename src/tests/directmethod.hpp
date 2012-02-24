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
dmtest (Connector<T>* rc) {

	Matrix<cxfl>  b1;  
	Matrix<cxfl>  tmxy;
	Matrix<float> tmz;  
	Matrix<float> r;   
	Matrix<float> b0;  
	
	Matrix<cxfl>  smxy;
	Matrix<float> smz;
	Matrix<float> roi;
	
	Matrix<float> g;   
	Matrix<float> j;
	
	std::string cf  = std::string (base + std::string(config));
	std::string odf = std::string (base + std::string("/simout.mat"));
	
	rc->ReadConfig (cf.c_str());
	
	std::string gf = std::string(base + std::string(rc->Attribute("gf"))); // gradient trajectories
	std::string pf = std::string(base + std::string(rc->Attribute("pf"))); // patterns
	std::string mf = std::string(base + std::string(rc->Attribute("mf"))); // maps
	
	// Gradients
#ifdef HAVE_MAT_H

	MXRead    (g, gf, rc->Attribute("g"));
	MXRead    (j, gf, rc->Attribute("j"));
	
	// Target excitation, ROI, sample
	MXRead    (r, pf, "r");
	MXRead (tmxy, pf, rc->Attribute("p"));
	tmz  = Matrix<float>::Zeros (r.Dim(1), 1);
	smxy = Matrix<cxfl>::Zeros  (r.Dim(1), 1);
	MXRead ( smz, pf, rc->Attribute("s"));
	MXRead ( roi, pf, rc->Attribute("roi"));

	// Maps
	MXRead ( b1, mf, rc->Attribute("b1"));
	MXRead ( b0, mf, rc->Attribute("b0"));
#endif	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}

	// Outgoing -------------
	
	rc->SetMatrix  (  "b1", b1  );
	rc->SetMatrix (   "g", g  );
	rc->SetMatrix (   "r", r   );
	rc->SetMatrix (  "b0", b0  );
	rc->SetMatrix  ("tmxy", tmxy);
	rc->SetMatrix ( "tmz", tmz );
	rc->SetMatrix  ("smxy", smxy);
	rc->SetMatrix ( "smz", smz );
	rc->SetMatrix (   "j", j   );
	rc->SetMatrix ( "roi", roi );
	// ---------------------
	
	rc->Process (test);
	
	// Incoming -------------
	
	Matrix<cxfl>  rf;
	Matrix<cxfl>  mxy;
	Matrix<float> mz;   

	rc->GetMatrix  ( "mxy", mxy);
	rc->GetMatrix (  "mz", mz);	
	rc->GetMatrix  ("tmxy", tmxy);
	rc->GetMatrix ( "tmz", tmz);	
	rc->GetMatrix  (  "rf", rf);

	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* od = matOpen (odf.c_str(), "w");

	if (od == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}

	MXDump  (mxy, od, "mxy");
	MXDump   (mz, od, "mz");
	MXDump (tmxy, od, "tmxy");
	MXDump  (tmz, od, "tmz");
	MXDump   (rf, od, "rf");

	if (matClose(od) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif

	return true;
	
}

