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

#include "CompressedSensing.hpp"
#include "Toolbox.hpp"
#include "DFT.hpp"

using namespace RRStrategy;


error_code
CompressedSensing::Init () {

	printf ("Intialising CompressedSensing ...\n");

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	int wli = 0;
	int m_fft = 0;

	Attribute ("tvw",     &m_csparam.tvw);
	Attribute ("xfmw",    &m_csparam.xfmw);
	Attribute ("l1",      &m_csparam.l1);
	Attribute ("pnorm",   &m_csparam.pnorm);
    printf ("  Weights: TV(%.2e) XF(%.2e) L1(%.2e)\n", m_csparam.tvw, m_csparam.xfmw, m_csparam.l1);
    printf ("  Pnorm: %.2e\n", m_csparam.pnorm);
	
	Attribute ("fft",     &m_csparam.fft);
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

	Attribute ("cgconv", &m_csparam.cgconv);
	Attribute ("cgiter", &m_csparam.cgiter);
	Attribute ("lsiter", &m_csparam.lsiter);
	Attribute ("lsa",    &m_csparam.lsa);
	Attribute ("lsb",    &m_csparam.lsb);
    printf ("  Iterations: CS(%i) CG(%i) LS(%i)\n", m_csiter, m_csparam.cgiter, m_csparam.lsiter);
	printf ("  Conv: CG(%.4f)\n", m_csparam.cgconv);
	printf ("  LS brackets: lsa(%.2e) lsb(%.2e)", m_csparam.lsa, m_csparam.lsb);

	

	m_initialised = true;
	printf ("... done.\n\n");

	return OK;

}


error_code
CompressedSensing::Process () {

	printf ("Processing CompressedSensing ...\n");
	ticks tic; 
	float ma;

	Matrix<cxfl>&  data  = Get<cxfl>   ("data");
	Matrix<float>& pdf   = Get<float>  ("pdf" );
	Matrix<float>& mask  = Get<float>  ("mask");
	Matrix<cxfl>&  pc    = Get<cxfl>   ("pc");
	Matrix<cxfl>&  im_dc = AddMatrix ("im_dc", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>  (data.DimVector())));
	Matrix<cxfl>   orig;

	printf ("  Geometry: %zuD (%zu,%zu,%zu)\n", ndims (data), 
		size(data,0), size(data,1), size(data,2));

	m_csparam.dwt = new DWT <cxfl> (data.Height(), wlfamily(m_wf), m_wm);

	/** -----  Which Fourier transform? **/
	m_csparam.ft  = (FT<float>*) new DFT<float> (size(data), mask, pc);

	
	// m_csparam.ft = (FT<float>*) new NFFT<float> (ms, M * shots, m, alpha);
	// m_csparam.ft = (FT<float>*) new NCSENSE<float> (sens, nk, m_cseps, m_csmaxit, m_lambda, m_fteps, m_ftmaxit);
	/*************************************/

	m_csparam.tvt = new TVOP ();

	FT<float>& dft = *m_csparam.ft;
	DWT<cxfl>& dwt = *m_csparam.dwt;
	
	im_dc    = data;
	im_dc   /= pdf;
	
	im_dc    = dft ->* im_dc;
	orig     = im_dc;
	
	ma       = max(abs(im_dc));
	im_dc   /= ma;
	data    /= ma;
	
	im_dc    = dwt * im_dc;
	
	printf ("  Running %i NLCG iterations ... \n", m_csiter); fflush(stdout);

	tic      = getticks();

	for (size_t i = 0; i < m_csiter; i++)
		NLCG (im_dc, data, m_csparam);

	printf ("  done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	im_dc    = dwt ->* im_dc * ma;
	data     = orig;

	return OK;

}


CompressedSensing::CompressedSensing() :
	m_wm(0), m_csiter(0), m_wf(0), m_dim(0) {}


CompressedSensing::~CompressedSensing() {}


error_code
CompressedSensing::Finalise() {return OK;}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CompressedSensing;
}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


