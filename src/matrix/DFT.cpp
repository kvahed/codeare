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
#include "Algos.hpp"
#include "CX.hpp"
#include "IO.hpp"

#include <fftw3.h>
#include <stdio.h>
#include <string.h>

template<> template<>
DFT<cxfl>::DFT (const size_t rank, const size_t sl, const Matrix<float> mask, const Matrix<cxfl> pc) :
	m_have_mask (false),
	m_have_pc (false) {
	
	int n[rank];
	
	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);

	MXDump (m_pc, "m_pc.mat", "m_pc");
	MXDump (m_cpc, "m_cpc.mat", "m_cpc");


	for (size_t i = 0; i < rank; i++)
		n[i]  = sl;

	m_N   = pow (sl, rank);

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}

template<> template<>
DFT<cxfl>::DFT (const Matrix<size_t>& size, const Matrix<float> mask, const Matrix<cxfl> pc) : m_N(1),
																						m_have_mask (false),
																						m_have_pc (false) {

	int rank = size.Size();
	int n[rank];

	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (size_t i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_ESTIMATE);

	m_initialised = true;

}


template<> 
DFT<cxfl>::DFT (const Matrix<size_t>& size) : m_N(1),
											 m_have_mask (false),
											 m_have_pc (false) {

	
	int rank = size.Size();
	int n[rank];

	for (size_t i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);
	m_out = (void*) fftwf_malloc (sizeof(fftwf_complex) * m_N);

	m_fwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplanf = fftwf_plan_dft (rank, n, (fftwf_complex*)m_in, (fftwf_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}


template<> template<>
DFT<cxdb>::DFT (const size_t rank, const size_t sl, const Matrix<double> mask, const Matrix<cxdb> pc) :
	m_have_mask (false),
	m_have_pc (false) {
	
	int n[rank];
	
	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (size_t i = 0; i < rank; i++)
		n[i]  = sl;

	m_N   = pow (sl, rank);

	m_in  = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);
	m_out = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);

	m_fwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_BACKWARD, FFTW_MEASURE);

	m_initialised = true;

}

template<> template<>
DFT<cxdb>::DFT (const Matrix<size_t>& size, const Matrix<double> mask, const Matrix<cxdb> pc) : m_N(1),
																						m_have_mask (false),
																						m_have_pc (false) {

	int rank = size.Size();
	int n[rank];

	if (mask.Size() > 1)
		m_have_mask = true;
	
	m_mask = mask;
	
	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);
	
	for (size_t i = 0; i < rank; i++) {
		n[i]  = (int)size[rank-1-i];
		m_N  *= n[i];
	}

	m_in  = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);
	m_out = (void*) fftw_malloc (sizeof(fftw_complex) * m_N);

	m_fwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_FORWARD,  FFTW_MEASURE);
	m_bwdplan = fftw_plan_dft (rank, n, (fftw_complex*)m_in, (fftw_complex*)m_out, FFTW_BACKWARD, FFTW_ESTIMATE);

	m_initialised = true;

}


template<>
DFT<cxfl>::~DFT () {

	fftwf_destroy_plan (m_fwdplanf);
	fftwf_destroy_plan (m_bwdplanf);

	fftwf_cleanup ();

	fftwf_free(m_in); 
	fftwf_free(m_out);

}


template<>
DFT<cxdb>::~DFT () {

	fftw_destroy_plan (m_fwdplan);
	fftw_destroy_plan (m_bwdplan);

	fftw_cleanup ();

	fftw_free(m_in); 
	fftw_free(m_out);

}


template<> Matrix<cxfl> 
DFT<cxfl>::Trafo (const Matrix<cxfl>& m) const {
	
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
DFT<cxfl>::Adjoint (const Matrix<cxfl>& m) const {

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
DFT<cxdb>::Trafo (const Matrix<cxdb>& m) const {
	
    Matrix<cxdb> res = fftshift(m);
	memcpy (m_in, &res[0], sizeof(fftw_complex) * m_N);
	if (m_have_pc)
		res *= m_pc;

	fftw_execute(m_fwdplan);

	memcpy (&res[0], m_out, sizeof(fftw_complex) * m_N);
	if (m_have_mask)
		res *= m_mask;

    return fftshift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxdb>
DFT<cxdb>::Adjoint (const Matrix<cxdb>& m) const {

    Matrix<cxdb> res = fftshift(m);
	if (m_have_mask)
		res *= m_mask;
	memcpy (m_in, &res[0], sizeof(fftw_complex) * m_N);

	fftw_execute(m_bwdplan);

	memcpy (&res[0], m_out, sizeof(fftw_complex) * m_N);
	if (m_have_pc)
		res *= m_cpc;

	return fftshift(res/sqrt((float)m.Size()));
	
}



