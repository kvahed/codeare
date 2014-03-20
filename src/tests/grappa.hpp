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


bool grappatest (Connector& rc) {

    using namespace codeare::matrix::io;
	// Incoming
	Matrix<cxdb> scan;   // Folded images O (Nc, RO, PE/AF, PE2)
	Matrix<cxdb> acs;   // Sensitivity maps O (Nc, RO, PE, PE2)
	
    IOContext ic (rc.GetElement("/config/data-in"), base_dir, READ);
    scan = ic.Read<cxdb>(rc.GetElement("/config/data-in/scan"));
    acs = ic.Read<cxdb>(rc.GetElement("/config/data-in/acs"));
    ic.Close();

	// Outgoing
	Matrix<cxdb> recon;   // Reconstructed unaliased O (RO, PE, PE2)

	if (rc.Init (test) != codeare::OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}
	
	// Prepare -------------
	
	rc.SetMatrix ("scan", scan); // Sensitivities
	rc.SetMatrix ("acs", acs); // Weights
	
	rc.Prepare   (test);
	
	// Process -------------
	
	rc.Process   (test);
	
	// Receive -------------
	
	//rc.GetMatrix ("recon", recon);  // Images
	
	// ---------------------
	
	rc.Finalise   (test);
	
	IOContext oc (rc.GetElement("/config/data-out"), base_dir, WRITE);
    oc.Write(recon, "recon");
    oc.Close();

	return true;

}

