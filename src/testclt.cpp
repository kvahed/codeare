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

#include "testclt.hpp"

int main (int argc, char** argv) {
	
	int  rc  = 0;

	if (init (argc, argv)) {
		
#ifdef LOCAL
		Connector<LocalConnector>*  con = new Connector<LocalConnector>  (name, verbose);
#else 
		Connector<RemoteConnector>* con = new Connector<RemoteConnector> (name, verbose);
#endif	
		
		std::string cf = std::string (base + std::string(config));
		con->ReadConfig (cf.c_str()); 

		if      (!strcmp (test, "CGSENSE"))               cgsensetest  (con);
		else if (!strcmp (test, "DirectMethod"))          dmtest       (con);
		else if (!strcmp (test, "SENSE"))                 sensetest    (con);
		else if (!strcmp (test, "NuFFT"))                 nuffttest    (con); 
		else if (!strcmp (test, "NuFFT_OMP"))             nuffttest    (con);
		else if (!strcmp (test, "GRAPPA"))                grappatest   (con);
		else if (!strcmp (test, "KTPoints"))              ktptest      (con);
		else if (!strcmp (test, "CompressedSensing"))     cstest       (con);
		else if (!strcmp (test, "mxtest"))                mxtest       (con);
		else if (!strcmp (test, "nitest"))                nitest       (con);
		else if (!strcmp (test, "iotest"))                iotest       (con);
		else if (!strcmp (test, "fftwtest"))              fftwtest     (con);
		else if (!strcmp (test, "dwttest"))               dwttest      (con);
		else if (!strcmp (test, "algotest"))              algotest     (con);
		else if (!strcmp (test, "syngotest"))             syngotest    (con);
		else if (!strcmp (test, "RelativeSensitivities")) resetest     (con);
		else if (!strcmp (test, "VDSpiral"))              vdspiraltest (con);
		else if (!strcmp (test, "KArb"))                  karbtest     (con);
		else if (!strcmp (test, "Creators"))              creatorstest (con);
		//else if (!strcmp (test, "MPI"))                   mpitest      (con);
		else if (!strcmp (test, "SHA256"))                sha256test   (con);
		else if (!strcmp (test, "ISMRMRD"))               ismrmrdtest  (con);
		//else if (!strcmp (test, "ocl"))                   oclmatrixtest(con);
		else                                              internaltest (con);

		delete con;
		
	} else
		
		rc = 1;	
	
	return rc;

}



