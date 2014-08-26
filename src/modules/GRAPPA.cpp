/*
 *  codeare Copyright (C) 2010-2011 Kaveh Vahedipour
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

#include "GRAPPA.hpp"
#include "Print.hpp"

using namespace RRStrategy;

codeare::error_code
GRAPPA::Init () {

	printf ("Intialising GRAPPA ...\n");

	Attribute ("nthreads",  &m_nthreads);
	Attribute ("lambda", &m_lambda);
	m_kernel_size = RHSList<size_t>("kernel_size");
	m_acceleration_factors = RHSList<size_t>("acceleration_factors");

	printf ("... done.\n\n");

	return codeare::OK;

}


codeare::error_code
GRAPPA::Prepare     () { 
	printf ("Preparing %s ...\n", Name());
	p.Set("nthreads", m_nthreads);
	p.Set("lambda", m_lambda);
	p.Set("kernel_size", m_kernel_size);
	p.Set("ac_data", Get<cxfl>("ac_data"));       // Sensitivities
    m_ft = CGRAPPA<float>(p);
    AddMatrix<cxfl> ("full_data");
	printf ("... done.\n\n");
	return codeare::OK;
}


codeare::error_code
GRAPPA::Process     () { 

	Matrix<cxfl>& undersampled_data = Get<cxfl>("undersampled_data");
	Matrix<cxfl>& full_data = Get<cxfl>("full_data");

	full_data = m_ft ->* undersampled_data;

	return codeare::OK;

}


codeare::error_code
GRAPPA::Finalise () {

	return codeare::OK;

}



// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new GRAPPA;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}
