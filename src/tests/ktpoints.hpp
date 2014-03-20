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

bool ktptest (Connector& con) {

    using namespace codeare::matrix::io;

	Matrix<cxfl>   target;
	Matrix<cxfl>   b1;
	Matrix<cxfl>   final;
	Matrix<float>  r;
	Matrix<float>  k;
	Matrix<float>  b0;
	
    IOContext ic (con.GetElement("/config/data-in"), base_dir, READ);
    target  = ic.Read<cxfl>(con.GetElement("/config/data-in/target"));
    b0      = ic.Read<float>(con.GetElement("/config/data-in/b0"));
    k	    = ic.Read<float>(con.GetElement("/config/data-in/k"));
    r	    = ic.Read<float>(con.GetElement("/config/data-in/r"));
    b1      = ic.Read<cxfl>(con.GetElement("/config/data-in/b1"));
    ic.Close();

	con.Init(test);

	con.SetMatrix ("target", target);
	con.SetMatrix ("b1",     b1);
	con.SetMatrix ("r",      r);
	con.SetMatrix ("k",      k);
	con.SetMatrix ("b0",     b0);
	
	con.Process    (test);

	Matrix<cxfl>   rf;
	Matrix<float>  grad;
	Matrix<float>  nrmse;

	con.GetMatrix ("nrmse",  nrmse);
	con.GetMatrix ("rf",     rf);
	con.GetMatrix ("grad",   grad);
	con.GetMatrix ("final", final);

    IOContext oc (con.GetElement("/config/data-out"), base_dir, WRITE);
    oc.Write (final, "final");
    oc.Write (rf,    "rf");
    oc.Write (grad,  "grad");
	oc.Write (nrmse, "nrmse");
	oc.Close();

	con.Finalise(test);

	return true;

}

