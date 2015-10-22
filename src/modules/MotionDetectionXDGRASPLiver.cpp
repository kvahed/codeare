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
#include "Lapack.hpp"
#include "Statistics.hpp"
#include "Smooth.hpp"
#include "LocalMaxima.hpp"

typedef TUPLE< Matrix<float>, Matrix<cxfl>, Matrix<float> > eig_t;

using namespace RRStrategy;

codeare::error_code MotionDetectionXDGRASPLiver::Init() {

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Prepare() {

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Process     () {

	Matrix<cxfl>& meas = Get<cxfl>("meas"), motion_signal_fft;
	Matrix<float> zip, tmp, si, cv, pc, v, motion_signal, motion_signal_new;
	Vector<float> f_x, res_peak, tmp_peak, res_peak_nor;
	Vector<size_t> idx, tmp_idx, fr_idx;
	float f_s, lf, hf;
	size_t nn, span = 5, pc_sel = 5;
	eig_t et;

	_nx = size(meas,0);
	_nv = size(meas,1);
	_nz = size(meas,2);
	_nc = size(meas,3);

	_ta = 95; // from raw data
	_tr = _ta/_nv;

	_time = _tr*linspace<float>(1,_nv,1);
	// Frequency stamp (only for the delay enhanced part)
	f_s = 1./_tr;
	f_x = linspace<float>(0,f_s,_nv).Container();
	f_x = f_x - .5*f_s; // frequency after FFT of the motion signal
	if (_nv/2%2==0)
	    f_x += f_x[_nv/4];

	nn  = 400; // Interpolation along z dimension
	// Take the central k-space points (c++ indexing)
	meas = squeeze(meas(CR(_nx/2),CR(),CR(),CR()));
	meas = zpad(meas,size(meas,0),nn,size(meas,2));
	meas = permute (meas,1,0,2);

	zip = flipud(abs(fftshift(fft(meas,0),0)));
	meas = permute (meas,1,0,2);
	// Remove some edge slices
	zip = zip(CR(21,size(zip,0)-40),CR(),CR());

	// Normalization the projection profiles
	for (size_t i = 0; i < _nc; ++i)
	    zip(R(),R(),R(i)) /= squeeze(repmat(mean(zip(CR(),CR(),CR(i)),0),size(zip,0),1));

	// Do PCA or SVD in each coil element to extract motion signal
	si  = permute (zip, 0, 2, 1);
	si  = transpose(resize(si, size(si,0)*_nc, _nv));
	cv  = cov(si);
	et  = eig2(cv);
	pc  = GET<0>(et);
	v   = abs(GET<1>(et));
	v   = v(CR(idx));
	pc  = pc(CR(),CR(idx));
	motion_signal = transpose(gemm(pc, si, 'C', 'C'));

	motion_signal_new = Matrix<float>(size(motion_signal,0),pc_sel);
	motion_signal_fft = motion_signal_new;
	for (size_t i = 0; i < pc_sel; ++i) {
		motion_signal_new(R(),R(i)) = smooth<float>(motion_signal(CR(),CR(i)),span); // TODO: smooth
		tmp = abs(fftshift(fft(motion_signal(CR(_nv/2+1,size(motion_signal,0)),CR(i))))); // TODO: fft (view)
		motion_signal_fft(R(),R(i)) = tmp/max(tmp(CR()));
	}
	// Take the component with the highest peak in respiratory motion range
	lf = 0.1; hf = 0.5; //Respiratory frequency range
	//tmp_idx = find(f_x>hf);
	//ft_idx=find(f_x<hf & f_x>lf);
/*
	tmp_peak = squeeze(motion_signal_fft(CR(tmp_idx),CR(),CR()));
	res_peak = squeeze(motion_signal_fft(CR(fr_idx),CR(),CR()));

	for (size_t i = 0; i < pc; ++i)
		res_peak_nor(R(),R(i)) = res_peak(R(),R(i))/max(tmp_peak(CR(),CR(i)));

	tt = max(res_peak_nor);
	tt = find(tt==max(tt));

	//Find the peak points
	t = 10;
	peaks = findpeaks (res_signal, 'MINPEAKDISTANCE', t);

	// Here the contrast enhancement curve need to be found and demodulated.
	ft = fittype( 'smoothingspline' );
	opts = fitoptions( ft );
	opts.SmoothingParam = 0.015;

	res_signal=motionsignal_new(CR(),CR(t));

	curve_data = prepareCurveData(Peak_Index,double(Res_Signal(Peak_Index)));
	fit_result = fit( curve_data.xData, curvedata.yData, ft, opts );

	cfval = coeffvalues(fit_result);
	ftmax = feval(CR(ft,cfval(1)),CR(1,ntviews));

	res_signal = res_signal-ftmax;
*/
	Add ("zip", zip);
	Add ("motion_signal", motion_signal_new);
	return codeare::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new MotionDetectionXDGRASPLiver;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}

