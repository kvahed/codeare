#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

void inufft2d (char* filename,int N,int M,int iteration, int weight) {

  int j,k,l;                    /* some indices  */
  double real,imag,t;           /* to read the real and imag part of a complex number */
  nfft_plan my_plan;            /* plan for the two dimensional nfft  */
  solver_plan_complex my_iplan; /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[2],my_n[2];          /* to init the nfft */
  double epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft*/
  int m = 6;
  double alpha = 2.0;
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=ceil(N*alpha);
  my_N[1]=N; my_n[1]=ceil(N*alpha);
  nfft_init_guru(&my_plan, 2, my_N, M, my_n, m, PRE_PHI_HUT| PRE_PSI|
                         MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                         FFTW_INIT| FFT_OUT_OF_PLACE,
                         FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* precompute lin psi if set */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);

  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  solver_init_advanced_complex(&my_iplan,(mv_plan_complex*)&my_plan, infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fin=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
    {
        fscanf(fin,"%le ",&my_iplan.w[j]);
    }
    fclose(fin);
  }

  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        int j2= j-N/2;
        int k2= k-N/2;
        double r=sqrt(j2*j2+k2*k2);
        if(r>(double) N/2)
          my_iplan.w_hat[j*N+k]=0.0;
        else
          my_iplan.w_hat[j*N+k]=1.0;
      }
    }
  }

  /* open the input file */
  fin=fopen(filename,"r");

  /* read x,y,freal and fimag from the knots */
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,"%le %le %le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1],
    &real,&imag);
    my_iplan.y[j] = real + _Complex_I*imag;
  }

	fclose(fin);

  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  /* precompute full psi */
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k]=0.0;

  t=nfft_second();

  /* inverse trafo */
  solver_before_loop_complex(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
      break;
    //fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
    //l+1,iteration);
    solver_loop_one_step_complex(&my_iplan);
  }


  t=nfft_second()-t;
  /*#ifdef HAVE_MALLINFO
  fprintf(stderr,"time: %e seconds mem: %i \n",t,nfft_total_used_memory());
#else
  fprintf(stderr,"time: %e seconds mem: mallinfo not available\n",t);
  #endif*/


  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {
    fprintf(fout_real,"%le ", creal(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,"%le ", cimag(my_iplan.f_hat_iter[k]));
  }

  fclose(fout_real);
  fclose(fout_imag);

  /* finalize the infft */
  solver_finalize_complex(&my_iplan);

  /* finalize the nfft */
  nfft_finalize(&my_plan);
}

void nufft2d (char * file, int N, int M)
{
  int j,k;            /* some variables */
  double real, imag;
  nfft_plan my_plan;  /* plan for the two dimensional nfft  */
  FILE* fp;
  FILE* fk;
  FILE* fi;

  /* initialise my_plan */
  nfft_init_2d(&my_plan,N,N,M);

  fp=fopen("knots.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1]);
  }
  fclose(fp);

  fi=fopen("input_f.dat","r");
  fk=fopen(file,"w");

  for(j=0;j<N;j++)
  {
    for(k=0;k<N;k++) {
		fscanf(fi, "%le %le ", &real,&imag);
		my_plan.f_hat[(N*j+k)] = real + imag * _Complex_I;
    }
  }

  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  nfft_trafo(&my_plan);

  for(j=0;j<my_plan.M_total;j++)
  {
    fprintf(fk,"%le %le %le %le\n",my_plan.x[2*j+0],my_plan.x[2*j+1],creal(my_plan.f[j]),cimag(my_plan.f[j]));
  }
  fclose(fk);
  fclose(fi);

  nfft_finalize(&my_plan);
}

void inufft3d (char* filename,int N,int M,int Z,int iteration, int weight) {
	
	int j,k,z,l;                  /* some variables  */
	double real,imag;             /* to read the real and imag part of a complex number */
	nfft_plan my_plan;            /* plan for the two dimensional nfft  */
	solver_plan_complex my_iplan;          /* plan for the two dimensional infft */
	FILE* fin;                    /* input file                         */
	FILE* fout_real;              /* output file (real part) */
	FILE* fout_imag;              /* output file (imag part) */
	int my_N[3],my_n[3];          /* to init the nfft */
	double epsilon=0.0000003;     /* tmp to read the obsolent z from 700.acs
									 epsilon is a the break criterion for
									 the iteration */
	unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft */
	
	/* initialise my_plan, specific.
	   we don't precompute psi */
	my_N[0]=Z; my_n[0]=ceil(Z*1.2);
	my_N[1]=N; my_n[1]=ceil(N*1.2);
	my_N[2]=N; my_n[2]=ceil(N*1.2);
	nfft_init_guru(&my_plan, 3, my_N, M, my_n, 6,
				   PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
				   MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE,
				   FFTW_MEASURE| FFTW_DESTROY_INPUT);
	
	/* precompute lin psi */
	if(my_plan.nfft_flags & PRE_LIN_PSI)
		nfft_precompute_lin_psi(&my_plan);
	
	if (weight)
		infft_flags = infft_flags | PRECOMPUTE_WEIGHT;
	
	/* initialise my_iplan, advanced */
	solver_init_advanced_complex(&my_iplan,(mv_plan_complex*)(&my_plan), infft_flags );
	
	/* get the weights */
	if(my_iplan.flags & PRECOMPUTE_WEIGHT) {
		
		fin=fopen("weights.dat","r");
		for(j=0;j<M;j++)
			fscanf(fin,"%le ",&my_iplan.w[j]);
		
		fclose(fin);
		
	}
	
	/* get the damping factors */
	if(my_iplan.flags & PRECOMPUTE_DAMP) {
		
		for(j=0;j<N;j++){
			for(k=0;k<N;k++) {
				for(z=0;z<N;z++) {
					int j2= j-N/2;
					int k2= k-N/2;
					int z2= z-N/2;
					double r=sqrt(j2*j2+k2*k2+z2*z2);
					if(r>(double) N/2)
						my_iplan.w_hat[z*N*N+j*N+k]=0.0;
					else
						my_iplan.w_hat[z*N*N+j*N+k]=1.0;
				}
			}
		}
	}
	
	/* open the input file */
	fin=fopen(filename,"r");
	
	/* open the output files */
	fout_real=fopen("output_real.dat","w");
	fout_imag=fopen("output_imag.dat","w");
	
	/* read x,y,freal and fimag from the knots */
	for(j=0;j<M;j++) {
		
		fscanf(fin,"%le %le %le %le %le ",&my_plan.x[3*j+1],&my_plan.x[3*j+2], &my_plan.x[3*j+0], &real,&imag);
		my_iplan.y[j] = real + _Complex_I*imag;
	}
	
	/* precompute psi */
	if(my_plan.nfft_flags & PRE_PSI)
		nfft_precompute_psi(&my_plan);
	
	/* precompute full psi */
	if(my_plan.nfft_flags & PRE_FULL_PSI)
		nfft_precompute_full_psi(&my_plan);
	
	/* init some guess */
	for(k=0;k<my_plan.N_total;k++)
		my_iplan.f_hat_iter[k]=0.0;
	
	/* inverse trafo */
	solver_before_loop_complex(&my_iplan);
	for(l=0;l<iteration;l++) {
		
		/* break if dot_r_iter is smaller than epsilon*/
		if(my_iplan.dot_r_iter<epsilon)
			break;
		fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
				l+1,iteration);
		solver_loop_one_step_complex(&my_iplan);
	}
	
	for(l=0;l<Z;l++) {
		for(k=0;k<N*N;k++)
			{
				/* write every Layer in the files */
				fprintf(fout_real,"%le ",creal(my_iplan.f_hat_iter[ k+N*N*l ]));
				fprintf(fout_imag,"%le ",cimag(my_iplan.f_hat_iter[ k+N*N*l ]));
			}
		fprintf(fout_real,"\n");
		fprintf(fout_imag,"\n");
	}
	
	fclose(fout_real);
	fclose(fout_imag);
	
	solver_finalize_complex(&my_iplan);
	nfft_finalize(&my_plan);
}

void nufft3d (char * file, int N, int M, int Z) {
	int j,k,l;                /* some variables */
	double real;
	nfft_plan my_plan;        /* plan for the three dimensional nfft  */
	FILE* fp,*fk;
	int my_N[3],my_n[3];      /* to init the nfft */
	
	
	/* initialise my_plan */
	//nfft_init_3d(&my_plan,Z,N,N,M);
	my_N[0]=Z; my_n[0]=ceil(Z*1.2);
	my_N[1]=N; my_n[1]=ceil(N*1.2);
	my_N[2]=N; my_n[2]=ceil(N*1.2);
	nfft_init_guru(&my_plan, 3, my_N, M, my_n, 6,
				   PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
				   MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE,
				   FFTW_MEASURE| FFTW_DESTROY_INPUT);
	
	fp=fopen("knots.dat","r");
	
	for(j=0;j<M;j++)
		fscanf(fp,"%le %le %le",&my_plan.x[3*(j)+1],
			   &my_plan.x[3*(j)+2],&my_plan.x[3*(j)+0]);
	
	fclose(fp);
	
	fp=fopen("input_f.dat","r");
	fk=fopen(file,"w");
	
	for(l=0;l<Z;l++) {
		for(j=0;j<N;j++) {
			for(k=0;k<N;k++) {
				//fscanf(fp,"%le ",&my_plan.f_hat[(N*N*(Z-l)+N*j+k+N*N*Z/2)%(N*N*Z)][0]);
				fscanf(fp,"%le ",&real);
				my_plan.f_hat[(N*N*l+N*j+k)] = real;
			}
		}
	}
	
    if(my_plan.nfft_flags & PRE_PSI)
		nfft_precompute_psi(&my_plan);
	
    nfft_trafo(&my_plan);
	
	
    for(j=0;j<my_plan.M_total;j++)
		fprintf(fk,"%le %le %le %le %le\n",my_plan.x[3*j+1],
				my_plan.x[3*j+2],my_plan.x[3*j+0],creal(my_plan.f[j]),cimag(my_plan.f[j]));
	
	
	
	fclose(fk);
	fclose(fp);
	
	nfft_finalize(&my_plan);
}

