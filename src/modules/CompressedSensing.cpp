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

#include "CompressedSensing.hpp"
#include "Toolbox.hpp"
#include "TVOP.hpp"
#include "IO.hpp"

using namespace RRStrategy;


RRSModule::error_code
CompressedSensing::Init () {

	printf ("Intialising CompressedSensing ...\n");

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	int wli = 0;
	int m_fft = 0;

	Attribute ("tvw",     &m_cgparam.tvw);
	Attribute ("xfmw",    &m_cgparam.xfmw);
	Attribute ("l1",      &m_cgparam.l1);
	Attribute ("pnorm",   &m_cgparam.pnorm);
    printf ("  Weights: TV(%.4f) XF(%.4f) L1(%.4f)\n", m_cgparam.tvw, m_cgparam.xfmw, m_cgparam.l1);
    printf ("  Pnorm: %f\n", m_cgparam.pnorm);
	
	Attribute ("fft",     &m_cgparam.fft);
	printf ("  FFT class: ");
	switch (m_fft) 
		{
		case 0:  printf ("%s", "Cartesian"); break;
		case 1:  printf ("%s", "Non-Cartesian"); break;
		default: printf ("%s", "Cartesian"); m_fft = 0; break;
		}
	printf ("\n");

	Attribute ("csiter", &m_csiter);
	Attribute ("wl_family", &m_wf);
	Attribute ("wl_member", &m_wm);
	printf ("  DWT(%i,%i)", m_wf, m_wm);
	
	if (m_wf < -1 || m_wf > 5)
		m_wf = -1;

	Attribute ("cgconv", &m_cgparam.cgconv);
	Attribute ("cgiter", &m_cgparam.cgiter);
	Attribute ("lsiter", &m_cgparam.lsiter);
	Attribute ("lsa",    &m_cgparam.lsa);
	Attribute ("lsb",    &m_cgparam.lsb);
    printf ("  Iterations: CS(%i) CG(%i) LS(%i)\n", m_csiter, m_cgparam.cgiter, m_cgparam.lsiter);
	printf ("  Conv: CG(%.4f)\n", m_cgparam.cgconv);
	printf ("  LS brackets: lsa(%.4f) lsb(%.4f)", m_cgparam.lsa, m_cgparam.lsb);

	

	m_initialised = true;
	printf ("... done.\n\n");

	return RRSModule::OK;

}


RRSModule::error_code
CompressedSensing::Process () {

	printf ("Processing CompressedSensing ...\n");
	ticks tic; 
	cxfl  ma;

	Matrix<cxfl>&   data  = GetCXFL   ("data");
	Matrix<double>& pdf   = GetRLDB   ("pdf" );
	Matrix<double>& mask  = GetRLDB   ("mask");
	Matrix<cxfl>&   pc    = GetCXFL   ("pc");
	Matrix<cxfl>&   im_dc = AddMatrix ("im_dc", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>  (data.Dim())));
	Matrix<cxfl>    orig;

    printf ("  Geometry: %iD (%lu,%lu,%lu)\n", HDim(data)+1, 
			data.Dim(0), data.Dim(1), data.Dim(2));

	m_cgparam.dwt = new DWT (data.Height(), wlfamily(m_wf));
	m_cgparam.dft = new DFT (HDim(data)+1, data.Height(), mask, pc);
	
	im_dc  = data;
	im_dc /= pdf;
	
	im_dc = m_cgparam.dft->Adjoint(im_dc);
	orig  = Matrix<cxfl>(im_dc);
	
	ma     = im_dc.Maxabs();
	im_dc /= ma;
	data  /= ma;
	
	im_dc  = m_cgparam.dwt->Trafo (im_dc);
	
	printf ("  Running %i NLCG iterations ... \n", m_csiter); fflush(stdout);

	tic    = getticks();

	for (int i = 0; i < m_csiter; i++)
		NLCG (im_dc, data, mask, m_cgparam);

	printf ("  done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	im_dc  = m_cgparam.dwt->Adjoint (im_dc) * ma;
	data   = orig;

	return OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CompressedSensing;
}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


