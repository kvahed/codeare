#include "CGSENSE.h"

CGSENSE::CGSENSE () {

	// Copy incoming data to m_temp             */
	Matrix < raw > m_temp = (m_raw);

	// Reshape for outgoing format
	Attribute ("Nx",             &m_raw.Dim(COL));
	Attribute ("Ny",             &m_raw.Dim(LIN));
	Attribute ("Nz",             &m_raw.Dim(SLC));

	for (int i = 3; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;

	// Keep iterations [LOTS OF DATA]? Allocate data.
	Attribute ("verbose", (int*) &m_verbose);
	if (m_verbose)
		m_raw.Dim(SET) = m_iter;

	m_raw.Reset();
	
	int N [2] = {m_raw.Dim(COL), m_raw.Dim(LIN)};
	int Nk[2] = {m_helper.Dim(COL), m_helper.Dim(LIN)}; 
	
	m_nufft.init (2, N, Nk, &m_helper[0], 2, 2.4, 32, 0, 0);


}


/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *
 * @param  in           Original discretised sample O (Nx x Ny)
 * @param  sens         Sensitivity maps            O (Nx x Ny x Nc)
 * @param  traj         k-space trajectory          O (Nk x 1)
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, noncart::strategy* ncs, Matrix<raw>* out) {
	
	/*
	  FT = cell(1,Nb_coils);
	  data_out = zeros(Nb_samples,Nb_coils);
	  %pos = [6 3 2 1 4 7 8 9];
	  
	  
	  %% Computation of full density k-space for each coil
	  for ind_coil = 1:Nb_coils
	  %figure(8), subplot (3,3,pos(ind_coil)), imagesc(abs(data_in.*sensitivity(:,:,ind_coil))), colormap gray, axis image;
	  FT{ind_coil} = fftshift(fftn(data_in.*sensitivity(:,:,ind_coil)));
	  end
	  % pause;
	  
	  %% reverse-gridding of full density k-space on actual k-trajectory
	  k_max(1) = abs(min(k(:,1)));
	  k_max(2) = abs(min(k(:,2)));
	  %figure(3);
	  for ind_sample = 1:Nb_samples
	  for ind_coil = 1:Nb_coils 
	  data_out(ind_sample,ind_coil) = FT{ind_coil}(k(ind_sample,1)+k_max(1)+1,k(ind_sample,2)+k_max(2)+1);
	  %subplot(3,3,ind_coil), imagesc(abs(data_out(ind_sample,ind_coil)));
	  end
	  end
	*/

	int nc   = sm->Dim(CHA);
	int nk   = in->Dim(COL); 

	int nx   = sm->Dim(COL);
	int ny   = sm->Dim(LIN);
	int nz   = sm->Dim(SLC);

	// Full density k-spaces 
	Matrix<raw>* tmp = new Matrix<raw> (nx, ny, nz, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); 
	
	((noncart::nufft*)ncs)->forward_2d(in, out);

		

		//fftshift(fftn(data_in.*sensitivity(:,:,ind_coil)));

	// Reverse grid full density k-space on actual trajectory
	return OK;

}


/**
 * @brief               Compute right hand side (i.e. Multiply EH, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sens         Sensitivity maps                 O (Nx x Ny x Nc)
 * @param  traj         k-space trajectory               O (Nk x 1)
 * @param  out          Returned product                 O (Nx x Ny)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, noncart::strategy* ncs, Matrix<raw>* out) {

	int ncoils   = sm->Dim(CHA);
	int nsamples = in->Size(); 

	return OK;

}


RRSModule::error_code
CGSENSE::Process () {

	Matrix<raw> *p, *s, *k, *q;

	// CG iterations
	for (int iter = 0; iter < m_iter; iter++) {

		EH (p, s, k, &m_nufft, q);

		delete p;
		p = new Matrix<raw> (*(q));

		E  (p, s, k, &m_nufft, q);

		/*
		  delta = r(:)'*r(:)/(a(:)'*a(:));
		  q     = eh(e(p , sensitivity, k), sensitivity, k);
		  b     = b + r(:)'*r(:)/(p(:)'*q(:))*p;
		  r_new = r - r(:)'*r(:)/(p(:)'*q(:))*q;
		  p     = r_new + r_new(:)'*r_new(:)/(r(:)'*r(:))*p;
		  r     = r_new;
		*/
		
	}

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

