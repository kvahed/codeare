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
#include "Interpolate.hpp"
#include "LocalMaxima.hpp"

using namespace RRStrategy;

codeare::error_code MotionDetectionXDGRASPLiver::Init() {

	try {
		_margin_top = GetAttr<float>("margin_top");
	} catch (const TinyXMLQueryException&) {}

	try {
		_margin_bottom = GetAttr<float>("margin_bottom");
	} catch (const TinyXMLQueryException&) {}

	try {
		_span = GetAttr<size_t>("smoothing_span");
	} catch (const TinyXMLQueryException&) {}

	try {
		_min_dist = GetAttr<size_t>("min_peak_dist");
	} catch (const TinyXMLQueryException&) {}

	try {
		_pc_sel = GetAttr<size_t>("n_pc_comp");
	} catch (const TinyXMLQueryException&) {}

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Prepare() {

	return codeare::OK;
}

codeare::error_code MotionDetectionXDGRASPLiver::Process     () {

	typedef float real_t;
	typedef cxfl  complex_t;

	Matrix<cxfl> meas = Get<cxfl>("meas");
	Matrix<float> zip, tmp, si, cv, pc, v, motion_signal, motion_signal_new, motion_signal_fft,
		res_peak, tmp_peak, res_peak_nor, tt, res_signal, ftmax;
	Vector<float> f_x;
	Vector<size_t> idx, tmp_idx, fr_idx, ft_idx, peaks;
	float f_s, lf, hf, ta, tr;
	size_t nn;
	eig_t<float> et;
    ta = wspace.PGet<float>("TA");
    tr = wspace.PGet<float>("TR");

    meas = squeeze(meas);
	std::cout << "  Incoming: " << size(meas) << std::endl;
	_nx = size(meas,0);
	_nc = size(meas,1);
	_nv = size(meas,2);
	_nz = size(meas,3);
	// Take the central k-space points (c++ indexing)
    meas = 1.0e6*squeeze(meas(CR(_nx/2),CR(),CR(),CR()));
    meas = permute(meas,1,2,0);
	std::cout << "  Reduced to center column and permuted: " << size(meas) << std::endl;
    
	std::cout << "  Analyse channel motion data ..." << std::endl;

	// Frequency stamp (only for the delay enhanced part)
	f_s = 1.0/(_nz*tr);
	f_x = linspace<float>(0,f_s,_nv/2).Container();
	f_x = f_x - .5*f_s; // frequency after FFT of the motion signal
	if (_nv/2%2==0)
	    f_x += f_x[_nv/4];

	nn  = 400; // Interpolation along z dimension
	meas = zpad(meas,size(meas,0),nn,size(meas,2));
	meas = permute (meas,1,0,2);

	zip = flipud(abs(fftshift(fft(meas,0),0)));
	meas = permute (meas,1,0,2);
	// Remove some edge slices
	zip = zip(CR(_margin_top,size(zip,0)-1-_margin_bottom),CR(),CR());

	std::cout << "  Normalise projection profiles ..." << std::endl;
	// Normalization the projection profiles
	for (size_t i = 0; i < _nc; ++i)
	    zip(R(),R(),R(i)) /= squeeze(repmat(mean(zip(CR(),CR(),CR(i)),0),size(zip,0),1));

	std::cout << "  Perform coil-wise PCA ..." << std::endl;
	// Do PCA or SVD in each coil element to extract motion signal
	si = permute (zip, 0, 2, 1);
	si = transpose(resize(si, size(si,0)*_nc, _nv));
	cv = cov(si);
	et = eigs(cv);
	pc = et.lv;
	v  = real(et.ev);
	v  = flipud(v);
	pc = fliplr(pc);
	motion_signal = transpose(gemm(pc, si, 'C', 'C'));

	std::cout << "  Choose component with highest respiratory amplitude ..." << std::endl;
	motion_signal_new = Matrix<float>(_nv  ,_pc_sel);
	motion_signal_fft = Matrix<float>(_nv/2,_pc_sel);
	for (size_t i = 0; i < _pc_sel; ++i) {
		motion_signal_new(R(),R(i)) = smooth(motion_signal(CR(),CR(i)),_span);
		tmp = abs(fftshift(fft(motion_signal(CR(_nv/2,size(motion_signal,0)-1),CR(i)),0,false)));
		motion_signal_fft(R(),R(i)) = tmp/max(tmp.Container());
	}

	std::cout << "  Detect peaks ..." << std::endl;
	// Take the component with the highest peak in respiratory motion range
	lf = 0.08; hf = 0.5; //Respiratory frequency range
	tmp_idx = find(f_x>hf);
	fr_idx=find(f_x<hf & f_x>lf);
	tmp_peak = squeeze(motion_signal_fft(CR(tmp_idx),CR()));
	res_peak = squeeze(motion_signal_fft(CR(fr_idx),CR()));
	res_peak_nor = res_peak;
	for (size_t i = 0; i < _pc_sel; ++i)
		res_peak_nor(R(),R(i)) /= mmax(motion_signal_fft(CR(),CR(i)));
	tt = max(res_peak_nor);
	std::cout << "  Respiration is component " << sort(tt,DESCENDING)[0] << " ..." << std::endl;
	res_signal = motion_signal_new(CR(),CR(sort(tt,DESCENDING)[0]));
	// Find peaks
	peaks = findLocalMaxima(res_signal,_min_dist,mean(res_signal)[0]);
	Matrix<float> peaks_i(peaks.size()+2,1), peaks_v = peaks_i;
	std::copy(peaks.begin(),peaks.end(),&peaks_i[1]);
	peaks_i[peaks.size()+1] = _nv+1.e-6;
	peaks_v(R(1,peaks.size()),R(0)) = res_signal(CR(peaks),CR(0));
	peaks_v[0] = max(res_signal[0],peaks_v[1]);
	peaks_v[peaks.size()+1] = peaks_v[peaks.size()];

	// Interpolate peaks
	std::cout << "  CSPLINE fit peaks and subtract..." << std::endl;
	ftmax = interp1(peaks_i, peaks_v, linspace<float>(0.,_nv-1,_nv));
	res_signal = res_signal-ftmax;

	Add ("ftmax", ftmax);
	Add ("motion_signal", motion_signal);
	Add ("motion_signal_new", motion_signal_new);
	Add ("motion_signal_fft", motion_signal_fft);
	Add ("res_peak", res_peak);
	Add ("res_peak_nor", res_peak_nor);
	Add ("peaks_i", peaks_i);
	Add ("peaks_v", peaks_v);
	Add ("res_signal", res_signal);
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

