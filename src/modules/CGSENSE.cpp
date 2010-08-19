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

}


RRSModule::error_code 
E  (Matrix<raw>* data_in, Matrix<raw>* sensitivity, Matrix<raw>* k, Matrix<raw>* data_out) {
	
	/*
	  Nb_coils = size(sensitivity,3);

	  Nb_samples = size(k,1);
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
	int ncoils   = sensitivity->Dim(CHA);
	int nsamples = data_in->Size(); 

	// Full density k-spaces 
	Matrix<raw> FT;

	// Reverse grid full density k-space on actual trajectory
	return OK;

}

RRSModule::error_code
EH (Matrix<raw>* data_in, Matrix<raw>* sensitivity, Matrix<raw>* k, Matrix<raw>* data_out) {

	int ncoils   = sensitivity->Dim(CHA);
	int nsamples = data_in->Size(); 

	return OK;

}

RRSModule::error_code
CGSENSE::Process () {

	Matrix<raw> *p, *s, *k, *q;

	// CG iterations
	for (int iter = 0; iter < m_iter; iter++) {

		EH (p, s, k, q);

		delete p;
		p = new Matrix<raw> (*(q));

		E  (p, s, k, q);

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

