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
#include "DFT.hpp"

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
	Attribute ("ntres", &m_ntres);
	m_params["ntres"] = m_ntres;
	m_nlines = 84/m_ntres;
	m_params["nlines"] = m_nlines;

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
	m_params["threads"] = (Attribute ("nthreads", &i) == TIXML_SUCCESS) ? i : 1;
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

	// data from scanner
	const Matrix<cxfl>& b1 = Get<cxfl>("sensitivities");
	const Matrix<float>& k = Get<float>("kspace");
	const Matrix<float>& w = Get<float>("weights");

	// sensitivity maps
	m_params["sensitivities"] = Get<cxfl>("sensitivities");

	// setup ft oper
	m_ft = NCSENSE<float>(m_params);
	m_ft.KSpace (Get<float>("kspace"));
	m_ft.Weights (Get<float>("weights"));

	// release RAM
	Free ("weights");
	Free("kspace");

	// prepare output
	Matrix<cxfl> img;
	Add ("image", img);

	return codeare::OK;
}



codeare::error_code XDGRASP::Process () {
    Matrix<cxfl> kdata = Get<cxfl>("signals");
    const Matrix<float>& traj = Get<float>("kspace");
    const Matrix<float>& resp = Get<float>("resp");
    
	size_t RO = 0, VIEW = 1, Z = 2, C = 3;
	// Cut edges
	kdata = fftshift(fft(kdata,1),1)/sqrt(size(kdata,2));
    kdata = kdata(Range(), Range(1,size(kdata,1)-1), Range());
	Vector<size_t> n = size(kdata);


	m_nt = (size_t)std::floor((float)n[VIEW]/
			(float)(unsigned_cast(m_params["ntres"])*unsigned_cast(m_params["nlines"])));
	m_params["nt"] = m_nt;

	// sort the data into two dynamic dimensions
	// one for contrast enhancement and one for respiration
	for (size_t i = 0; i < m_nt; ++i) {
		Matrix<cxfl> kdata_under = kdata(Range(), Range(i*m_ntres*m_nlines,(i+1)*m_ntres*m_nlines), Range());
		Matrix<float> traj_under = traj (Range(), Range(i*m_ntres*m_nlines,(i+1)*m_ntres*m_nlines));
		Matrix<float> resp_under = resp (Range(         i*m_ntres*m_nlines,(i+1)*m_ntres*m_nlines));
	}
	for (size_t i = 0; i < m_nt; ++i) {

	}
	/*for ii=1:nt
	    kdata_Under(:,:,:,:,ii)=kdata(:,(ii-1)*ntres*nline+1:ii*ntres*nline,:,:);
	    Traj_Under(:,:,ii)=Traj(:,(ii-1)*ntres*nline+1:ii*ntres*nline);
	    DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*ntres*nline+1:ii*ntres*nline);
	    Res_Signal_Under(:,ii)=Res_Signal((ii-1)*ntres*nline+1:ii*ntres*nline);
	end
	for ii=1:nt
	    tmp1=kdata_Under(:,:,:,:,ii);
	    tmp2=Traj_Under(:,:,ii);
	    tmp3=DensityComp_Under(:,:,ii);
	    [~,index]=sort(Res_Signal_Under(:,ii),'descend');
	    tmp1=tmp1(:,index,:,:);
	    tmp2=tmp2(:,index,:);
	    tmp3=tmp3(:,index);
	    for jj=1:ntres
	        kdata_Under1(:,:,:,:,jj,ii)=tmp1(:,(jj-1)*nline+1:jj*nline,:,:);
	        Traj_Under1(:,:,jj,ii)=tmp2(:,(jj-1)*nline+1:jj*nline);
	        DensityComp_Under1(:,:,jj,ii)=tmp3(:,(jj-1)*nline+1:jj*nline);
	    end
	end

	param.y=squeeze(double(kdata_Under1));
	b1=squeeze(double(b1));b1=b1/max(abs(b1(:)));
	param.E=MCNUFFT3D_MP(squeeze(Traj_Under1),squeeze(DensityComp_Under1(:,:,:,:)),b1);
	recon_cs=param.E'*param.y;

	param.TV=TV_Temp4D; % TV constraint along contrast dimension
	param.W=TV_Temp3DRes; % % TV constraint along respiratory dimension
	param.TVWeight=max(abs(recon_cs(:)))*0.03;
	param.L1Weight=max(abs(recon_cs(:)))*0.01;
	param.nite = 8;param.display=1;


	for n=1:3
	    recon_cs = CSL1NlCg(recon_cs,param);
	end
	*/
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
