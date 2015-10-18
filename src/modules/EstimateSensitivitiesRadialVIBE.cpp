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

#include "Algos.hpp"
#include "EstimateSensitivitiesRadialVIBE.hpp"

using namespace RRStrategy;

codeare::error_code EstimateSensitivitiesRadialVIBE::Init () {

	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Prepare () {
	Matrix<cxfl> sensitivities;
	Add<cxfl>("sensitivities", sensitivities);
	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Process () {
	Matrix<cxfl>& meas = Get<cxfl>("meas");
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");
	meas = squeeze(meas);
	meas = meas(CR(),CR(),CR(),CR(1,size(meas,3)-2));

	Vector<size_t> perm(4);
	perm[0]=0; perm[1]=3; perm[2]=2; perm[3]=1;

	sensitivities = permute(meas,perm);
	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Finalise () {
	return codeare::OK;
}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new EstimateSensitivitiesRadialVIBE;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}
