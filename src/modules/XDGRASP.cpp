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

#include "XDGRASP.hpp"

using namespace RRStrategy;


codeare::error_code XDGRASP::Init () {

	size_t s;
	bool b;
	float f;
	int i;

	std::cout << "Intialising XDGRASP ..." << std::endl;

	// XD GRASP specific
	Attribute ("nx", &s);
	m_params["nx"] = s;
	Attribute ("ntres", &s);
	m_params["ntres"] = s;
	printf ("  Image side length (%lu) Resp phases (%lu)\n",
			unsigned_cast(m_params["nx"]), unsigned_cast(m_params["ntres"]));

	// NC SENSE specific
	assert (Attribute ("nk", &s) == TIXML_SUCCESS);
	m_params["nk"] = s;
	m_params["ftiter"]  = (Attribute ("ftmaxit", &s) == TIXML_SUCCESS) ? s : 3;
	m_params["verbose"] = (Attribute ("verbose", &i) == TIXML_SUCCESS) ? i : 0;
	m_params["cgiter"]  = (Attribute ("cgmaxit", &s) == TIXML_SUCCESS) ? s : 10;
	m_params["cgeps"]   = (Attribute ("cgeps",   &f) == TIXML_SUCCESS) ? f : 0.;
	m_params["lambda"]  = (Attribute ("lambda",  &f) == TIXML_SUCCESS) ? f : 0.;
	m_params["3rd_dim_cart"] = (Attribute ("cart_3rd_dim", &b) == TIXML_SUCCESS) ? b : false;
	m_params["m"]       = (Attribute ("m",       &s) == TIXML_SUCCESS) ? s : 1;

#pragma omp parallel default (shared)
	{
		if (omp_get_thread_num() == 0)
			m_params["np"] = (Attribute ("nthreads", &s) == TIXML_SUCCESS) ? s : omp_get_num_threads();
	}

	std::cout << "... done." << std::endl;
	return codeare::OK;
}


codeare::error_code XDGRASP::Prepare () {
	const Matrix<cxfl>& b1 = Get<cxfl>("sensitivities");
	const Matrix<float>& k = Get<float>("kspace");
	const Matrix<float>& w = Get<float>("weights");

	m_params["sensitivities"] = Get<cxfl>("sensitivities");

	m_ft = NCSENSE<float>(m_params);

	m_ft.KSpace (Get<float>("kspace"));
	m_ft.Weights (Get<float>("weights"));
	Free ("weights");
	Free("kspace");

	return codeare::OK;
}


codeare::error_code XDGRASP::Process () {
	const Matrix<cxfl>& kdata = Get<cxfl>("signals");
	//kdata = permute (kdata());
	return codeare::OK;
}


codeare::error_code XDGRASP::Finalise () {
	return codeare::OK;
}


// the class factory
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new XDGRASP;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}
