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

#include "Creators.hpp"
#include "MotionDetectionXDGRASPLiver.hpp"
#include "DFT.hpp"

using namespace RRStrategy;

codeare::error_code MotionDetectionXDGRASPLiver::Init() {

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Prepare() {

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Process     () {

	Matrix<cxfl> kdata = Get<cxfl>("kdata"), zip;
	Vector<float> f_x;
	float f_s;
	size_t nn;

	_nx = size(kdata,0);
	_ntviews = size(kdata,1);
	_nz = size(kdata,2);
	_nc = size(kdata,3);

	_ta = 95; // from raw data
	_tr = _ta/_ntviews;

	_time = (linspace<float>(1,_ntviews,1)).Container()*_tr;

	// Frequency stamp (only for the delay enhanced part)
	f_s = 1./_tr;
	f_x = linspace<float>(0,f_s,_ntviews).Container();
	f_x = f_x - .5*f_s; // frequency after FFT of the motion signal
	if (_ntviews/2%2==0)
	    f_x += f_x[_ntviews/4];

	nn = 400; // Interpolation along z dimension
	// Take the central k-space points
	kdata = squeeze(kdata(CR(_nx/2+1),CR(),CR(),CR()));
	kdata = zpad(kdata,size(kdata,0),400,size(kdata,2));
	zip   = abs(fftshift(fft(squeeze(kdata),2),2));
	zip   = flipud(zip);

	// Remove some edge slices
	zip   = zip(CR(21,size(zip,0)-40),CR(),CR());

	Add ("zip", zip);

	return codeare::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new MotionDetectionXDGRASPLiver;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}

