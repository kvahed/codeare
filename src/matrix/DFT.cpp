/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "DFT.hpp"
#include "FFT.hpp"
#include "Algos.hpp"
#include "CX.hpp"
#include "IO.hpp"

#include <fftw3.h>
#include <stdio.h>
#include <string.h>


DFT::DFT (const size_t rank, const size_t sl, const Matrix<double> mask, const Matrix<cxfl> pc) :
	m_have_mask (false),
	m_have_pc (false) {
	
	int n[rank];
	
	if (mask.Size() > 1) {
		m_have_mask = true;
		m_mask = mask;
	}
	
	if (pc.Size() > 1) {
		m_have_pc = true;
		m_pc      = pc;
		m_cpc     = conj(pc);

		MXDump (m_pc, "m_pc.mat", "m_pc");
		MXDump (m_cpc, "m_cpc.mat", "m_cpc");
	}

	for (size_t i = 0; i < rank; i++)
		n[i]  = sl;

	m_N   = pow (sl, rank);

	m_in  = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, m_in, m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, m_in, m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}


DFT::DFT (const Matrix<int> size, const Matrix<double> mask, const Matrix<cxfl> pc) : m_N(1),
																						m_have_mask (false),
																						m_have_pc (false) {
	
	int rank = size.Size();
	int n[rank];

	
	if (mask.Size() > 1) {
		m_have_mask = true;
		m_mask = mask;
	}
	
	if (pc.Size() > 1) {
		m_have_pc = true;
		m_pc = pc;
		m_cpc = conj(pc);
	}

	for (size_t i = 0; i < rank; i++) {
		n[i]  = size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, m_in, m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, m_in, m_out, FFTW_BACKWARD, FFTW_ESTIMATE);

	m_initialised = true;

}


DFT::~DFT () {

	fftwf_destroy_plan (m_fwdplanf);
	fftwf_destroy_plan (m_bwdplanf);

	fftwf_cleanup ();

	fftwf_free(m_in); 
	fftwf_free(m_out);

}


template<> Matrix<cxfl> 
DFT::Trafo (const Matrix<cxfl>& m) const {
	
    Matrix<cxfl> res = fftshift(m);
	memcpy (m_in, &res[0], sizeof(fftwf_complex) * m_N);
	if (m_have_pc)
		res *= m_pc;

	fftwf_execute(m_fwdplanf);

	memcpy (&res[0], m_out, sizeof(fftwf_complex) * m_N);
	if (m_have_mask)
		res *= m_mask;

    return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxfl>
DFT::Adjoint (const Matrix<cxfl>& m) const {

    Matrix<cxfl> res = fftshift(m);
	if (m_have_mask)
		res *= m_mask;
	memcpy (m_in, &res[0], sizeof(fftwf_complex) * m_N);

	fftwf_execute(m_bwdplanf);

	memcpy (&res[0], m_out, sizeof(fftwf_complex) * m_N);
	if (m_have_pc)
		res *= m_cpc;

	return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxdb> 
DFT::Trafo (const Matrix<cxdb>& m) const {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p;

	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

    return fftshift(res / (float)res.Size());
	
}


template<> Matrix<cxdb>
DFT::Adjoint (const Matrix<cxdb>& m) const {

	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p; 
	
	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

	return fftshift(res);
	
}


template<> Matrix<cxfl> 
DFT::operator* (const Matrix<cxfl>& m) const {

	return Trafo(m);

}


template<> Matrix<cxfl> 
DFT::operator->* (const Matrix<cxfl>& m) const {

	return Adjoint (m);

}


