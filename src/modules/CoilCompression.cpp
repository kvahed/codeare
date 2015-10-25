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

#include "CoilCompression.hpp"
#include "Algos.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

using namespace RRStrategy;
static const char ECON = 'S';

codeare::error_code CoilCompression::Init () {
	_coil_dimension = GetAttr<size_t>("coil_dimension");
	_coils_left = GetAttr<size_t>("coils_remaining");
	return codeare::OK;
}

codeare::error_code CoilCompression::Prepare () {

	return codeare::OK;
}

codeare::error_code CoilCompression::Process () {

	typedef TUPLE<Matrix<cxfl>,Matrix<float>,Matrix<cxfl> > svd_t;

	Matrix<cxfl>& meas = Get<cxfl> ("meas"), V;
	Matrix<float> S;
	meas = squeeze(meas);

	// Permute coils to outermost dimension
	Vector<size_t> dims = size(meas), order(dims.size());
	std::iota(order.begin(), order.end(), 0);
	order.erase(order.begin()+_coil_dimension);
	order.push_back(_coil_dimension);
	std::cout << "  Incoming: " << dims << std::endl;
	meas = permute(meas,order);
	dims = size(meas);
	size_t ncoils = dims.back();
	std::cout << "  Permuted: " << dims << std::endl;
	std::cout << "  #Coils: " << ncoils << std::endl;

	// SVD measurement data
	std::cout << "  Performing SVD ..." << std::endl;
	meas = resize(meas, numel(meas)/ncoils, ncoils);
	svd_t usv = svd2 (meas, ECON); V = GET<2>(usv);
	std::cout << "  Recombining virtual coils ..." << std::endl;

	// Recombine compressed data
	V = V(CR(),CR(0,_coils_left-1));
	meas = meas->*V;

	// Resize data
	dims.back() = _coils_left;
	meas = resize(meas,dims);
	std::cout << "  Outgoing: " << size(meas) << std::endl;

	return codeare::OK;
}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new CoilCompression;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}
