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

    int error =1;
	
	if (init (argc, argv)) {
		
#ifdef LOCAL
		Connector<LocalConnector>*  con = new Connector<LocalConnector>  (name, verbose);
#else 
		Connector<RemoteConnector>* con = new Connector<RemoteConnector> (name, verbose);
#endif	
		
		std::string cf = base;
		cf += "/"; 
		cf += std::string(config);
		if (cf.compare(std::string(base)) != 0)
            con->ReadConfig (cf.c_str());

		if      (!strcmp (test, "CGSENSE"))               cgsensetest  (con);
		else if (!strcmp (test, "SENSE"))                 sensetest    (con);
		else if (!strcmp (test, "NuFFT"))                 nuffttest    (con); 
		else if (!strcmp (test, "NuFFT_OMP"))             nuffttest    (con);
		else if (!strcmp (test, "GRAPPA"))                grappatest   (con);
		else if (!strcmp (test, "KTPoints"))              ktptest      (con);
		else if (!strcmp (test, "CompressedSensing"))     cstest       (con);
#ifndef _MSC_VER
        else if (!strcmp (test, "VDSpiral"))              vdspiraltest (con);
		else if (!strcmp (test, "KArb"))                  karbtest     (con);
        else if (!strcmp (test, "SliceGrad"))             slicegrad    (con);
#endif
        else if (!strcmp (test, "DummyRecon"))            dummytest    (con);
        /*
		else if (!strcmp (test, "DirectMethod"))          dmtest       (con);
		else if (!strcmp (test, "RelativeSensitivities")) resetest     (con);
        */
        else    {printf ("ERROR: Cannot call unknow module %s\n", test); return 1;}

        error = 0;

	} 
		
    return error;	

}



