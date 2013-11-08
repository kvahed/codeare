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
ktptest (Connector<T>* rc) {

    using namespace codeare::matrix::io;

	Matrix<cxfl>   target;
	Matrix<cxfl>   b1;
	Matrix<cxfl>   final;
	Matrix<float>  r;
	Matrix<float>  k;
	Matrix<float>  b0;
	
    IOContext ic (rc->GetElement("/config/data-in"), base, READ);
    target  = ic.Read<cxfl>(rc->GetElement("/config/data-in/target"));
    b0      = ic.Read<float>(rc->GetElement("/config/data-in/b0"));
    k	    = ic.Read<float>(rc->GetElement("/config/data-in/k"));
    r	    = ic.Read<float>(rc->GetElement("/config/data-in/r"));
    b1      = ic.Read<cxfl>(rc->GetElement("/config/data-in/b1"));
    ic.Close();

	rc->Init(test);

	rc->SetMatrix ("target", target);
	rc->SetMatrix ("b1",     b1);
	rc->SetMatrix ("r",      r);
	rc->SetMatrix ("k",      k);
	rc->SetMatrix ("b0",     b0);
	
	rc->Process    (test);

	Matrix<cxfl>   rf;
	Matrix<float>  grad;
	Matrix<float>  nrmse;

	rc->GetMatrix ("nrmse",  nrmse);
	rc->GetMatrix ("rf",     rf);
	rc->GetMatrix ("grad",   grad);
	rc->GetMatrix ("final", final);

	rc->Finalise(test);

	
    IOContext oc (rc->GetElement("/config/data-out"), base, WRITE);
    oc.Write (final, "final");
    oc.Write (rf,    "rf");
    oc.Write (grad,  "grad");
	oc.Write (nrmse, "nrmse");
	oc.Close();

	return true;

}

