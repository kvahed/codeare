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

#include "matrix/Creators.hpp"
#include "matrix/Algos.hpp"

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
	Matrix<float> t;   
	Matrix<float> j;
	
	std::string cf  = std::string (base + std::string(config));
	std::string odf = std::string (base + std::string("/simout.mat"));
	
	std::string b1_map_name;
	std::string target_name;
	std::string b0_map_name;
	std::string

    IOContext f (rc->GetElement("/config/data-in"), base, READ);
    g = fread (f, "")
	
#ifdef HAVE_MAT_H

	// Gradients
	Read    (g, rc->GetElement("/config/data/g"),    base);
	Read    (j, rc->GetElement("/config/data/j"),    base);
	Read    (t, rc->GetElement("/config/data/t"),    base);
	
	// Target excitation, ROI, sample
	Read    (r, rc->GetElement("/config/data/r"),    base);
	Read (tmxy, rc->GetElement("/config/data/tmxy"), base);
	tmz  = zeros<float> (size(r,1), 1);
	smxy = zeros<cxfl>  (size(r,1), 1);
	Read ( smz, rc->GetElement("/config/data/smz"), base);
	Read ( roi, rc->GetElement("/config/data/roi"), base);
	
	// Maps
	Read ( b1, rc->GetElement("/config/data/b1"),   base);
	Read ( b0, rc->GetElement("/config/data/b0"),   base);

#endif	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}

	// Outgoing -------------
	
	rc->SetMatrix  ( "b1", b1  );
	rc->SetMatrix (   "g", g  );
	rc->SetMatrix (   "t", t  );
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

	h5write (mxy, od "mxy");
	h5write (mz, od, "mz");
	h5write (tmxy, od, "tmxy");
	h5write (tmz, od, "tmz");
	h5write (rf, od, "rf");

	if (matClose(od) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif

	return true;
	
}

