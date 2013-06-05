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


# ifndef __DWT2_HPP__

# define __DWT2_HPP__


/**
 * @brief  Supported wavelet families
 */
//enum wlfamily {
//
//	ID = -1,                  /**< Identity transform*/
//	WL_DAUBECHIES,
//	WL_DAUBECHIES_CENTERED,
//	WL_HAAR,
//	WL_HAAR_CENTERED,
//	WL_BSPLINE,
//	WL_BSPLINE_CENTERED

//};


# include "Matrix.hpp"
# include "DWT.hpp"



/**
 * @brief 2D Discrete wavelet transform for Matrix template (from GSL)
 */
template <class T>
class DWT2 {


//	typedef typename DWTTraits<T>::Type Type;


public:


	/**
	 * @brief Construct 2D Wavelet transform with wavelet class and side length
	 *
	 * @param  sl      Side length
	 * @param  wf      Wavelet family (default none, i.e. ID)
	 * @param  wm      Familty member (default 4)
	 */
    DWT2 (const int dim = 2)
		: m_dim (dim),
//		  m_lpf_d (2),
//		  m_lpf_r (2),
//		  m_hpf_d (2),
//		  m_hpf_r (2)
		  m_lpf_d (8),    // (db4)
		  m_lpf_r (8),
		  m_hpf_d (8),     // (db4)
		  m_hpf_r (8)
	{

    	float norm_factor = 1 / sqrt (m_dim);

    	// set high pass filter (haar)
//    	m_lpf_d [0] = 1; m_lpf_d [1] = 1;
//    	m_lpf_d *= norm_factor;
//        m_lpf_r [0] = 1; m_lpf_r [1] = 1;
//        m_lpf_r *= norm_factor;
    	m_lpf_d [0] = -0.0106; m_lpf_d [1] = 0.0329; m_lpf_d [2] = 0.0308; m_lpf_d [3] = -0.1870;
    	m_lpf_d [4] = -0.0280; m_lpf_d [5] = 0.6309; m_lpf_d [6] = 0.7148; m_lpf_d [7] = 0.2304;
        m_lpf_r [0] = 0.2304; m_lpf_d [1] = 0.7148; m_lpf_d [2] = 0.6309; m_lpf_d [3] = -0.0280;
        m_lpf_r [4] = -0.1870; m_lpf_d [5] = 0.0308; m_lpf_d [6] = 0.0329; m_lpf_d [7] = -0.0106;


    	// set low pass filter (haar)
//    	m_hpf_d [0] = -1; m_hpf_d [1] = 1;
//    	m_hpf_d *= norm_factor;
//        m_hpf_r [0] = 1; m_hpf_r [1] = -1;
//        m_hpf_r *= norm_factor;
    	m_hpf_d [0] = -0.2304; m_hpf_d [1] = 0.7148; m_hpf_d [2] = -0.6409; m_hpf_d [3] = -0.0280;
    	m_hpf_d [4] = 0.1870; m_hpf_d [5] = 0.0308; m_hpf_d [6] = -0.0329; m_hpf_d [7] = -0.0106;
        m_hpf_d [0] = -0.0106; m_hpf_d [1] = -0.0329; m_hpf_d [2] = 0.0308; m_hpf_d [3] = 0.1870;
        m_hpf_d [4] = -0.0280; m_hpf_d [5] = -0.6309; m_hpf_d [6] = 0.7148; m_hpf_d [7] = -0.2304;


    }
    

    virtual
    ~DWT2 () {
        

        
    }
    
    
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline
    Matrix<T>
	Trafo        (const Matrix<T>& m) {
		return TransformForward (m);
	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline
	Matrix<T>
	Adjoint      (const Matrix<T>& m) {
		return TransformBackwards (m);
	}
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline
    Matrix<T>
	operator*    (const Matrix<T>& m) {
		return Trafo(m);
	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline
    Matrix<T>
	operator->* (const Matrix<T>& m) {
		return Adjoint(m);
	}
	
	
private:
	

    /**
     * variable definitions
     */

    // dimension of DWT2
    const int m_dim;

    // low pass filters
    Matrix <double> m_lpf_d;
    Matrix <double> m_lpf_r;

    // high pass filters
    Matrix <double> m_hpf_d;
    Matrix <double> m_hpf_r;


    /**
     * function definitions
     */

	/**
	 * @brief   Transform
	 *
	 * @param   m   To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	inline
	Matrix<T>
	TransformForward    	(const Matrix<T>& m)
	{
        
        Matrix<T> res (m.DimVector());

        if (m_dim == 2)
        {

        	// create temporary memory
        	T * temp = (T*) malloc(3*m.Height()*m.Width()*sizeof(T));

        	// create vars from mex function
        	int J = 0, nn;
        	for (nn = 1; nn < m.Height (); nn *= 2 )
        		J ++;
        	if (nn  !=  m.Height ()){
        		std::cout << "FWT2 requires dyadic length sides" << std::endl;
        		assert (false);
        	}
        	const int ell = 0;

        	// call dpwt2
        	res = dpwt2 (m, ell, J, temp);

        	free (temp);

        }

        return res;
        
    }


	/**
	 * backwards transform
	 */
	inline
	Matrix <T>
	TransformBackwards		(const Matrix <T> & m)
	{

		Matrix <T> res (m.DimVector());

        if (m_dim == 2)
        {

        	// create temporary memory
        	T * temp = (T*) malloc (4*m.Height()*m.Width()*sizeof(T));

        	// create vars from mex function
        	int J = 0, nn;
        	for (nn = 1; nn < m.Height (); nn *= 2 )
        		J ++;
        	if (nn  !=  m.Height ()){
        		std::cout << "IWT2 requires dyadic length sides" << std::endl;
        		assert (false);
        	}
        	const int ell = 0;

        	// call dpwt2
        	res = idpwt2 (m, ell, J, temp);

        	free (temp);

        }

		return res;

	}


	/**
	 *  COPIED FROM WAVELAB IMPLEMENTATION
	 */

	Matrix <T>
	dpwt2		(const Matrix <T> & sig, const int ell, const int J, T * temp)
	{

		Matrix <T> res (sig.DimVector ());

		T *wcplo,*wcphi,*templo,*temphi;
		int k,j,nj;
//		copydouble(&sig[0],&res[0],sig.Height()*sig.Width());

		// assign signal to result matrix
		res = sig;

		templo = &temp[sig.Height ()];
		temphi = &temp[2*sig.Height()];

		int nr = sig.Height();
		int nc = sig.Width();
		nj = sig.Height();
		for( j=(J-1); j >= ell; --j){
			for( k=0; k < nj; k++){
				wcplo = &res[k*nr];
				wcphi = &res[k*nr + nj/2];
				copydouble(wcplo,temp,nj);
				downlo(temp, nj,wcplo);
				downhi(temp, nj,wcphi);
			}
			for( k=0; k < nj; k++){
				unpackdouble(&res[0],nj,nc,k,temp);
				downlo(temp, nj,templo);
				downhi(temp, nj,temphi);
				packdouble(templo,nj/2,nc,k,&res[0]);
				packdouble(temphi,nj/2,nc,k,&res[nj/2*nr]);
			}
			nj = nj/2;
		}

		return res;

	}

	void
	unpackdouble	(const T * const x, const int n, const int nc, const int k, T * const y)
	{
		int i;
		for( i=0; i < n; i++)
			y[i] = x[k+nc*i];
	}

	void
	packdouble		(const T * const x, const int n, const int nc, const int k, T * const y)
	{
		int i;
		for( i=0; i < n; i++)
			y[k+nc*i] = x[i];
	}


	void
	copydouble		(const T * const x, T * const y, const int n)
	{

		// TODO: use operator= instead ...
		for (int i = n; i > 0; --i)
		{
			y[i] = x[i];
		}

	}


	void
	downhi			(const T * const x, const int n, T * const y)
	{
		int n2, mlo, i, h, j;
		T s;

		int m = m_hpf_d.Dim (0);

		/* highpass version */
		n2 = n/2;
		mlo = m/2-1;
		if(2*mlo+1 - (m-1) < 0) mlo++;
		for( i= mlo; i<n2; i++) {
			s = 0.;
			for( h=0; h < m; h++)
				s += ((T)m_hpf_d[h])*x[2*i+1-h];
			y[i] = s;
		}
		if(mlo > n2) mlo = n2;
		/* fix up edge values */
		for( i= 0; i<mlo; i++) {
			s = 0.;
			j = 2*i+1;
			for( h=0; h < m; h++) {
				if(j < 0) j += n;
				s += ((T)m_hpf_d[h])*x[j];
				--j;
			}
			y[i] = s;
		}
	}


	void
	downlo		(const T * const x, const int n, T * const y)
	{

		int n2, mlo, mhi, i, h, j;
		T s;

		int m = m_lpf_d.Dim (0);

		/*lowpass version */
		n2 = n/2;
		mlo = m /2;
		mhi = n2-mlo;
		if(2*mhi + (m-1) >= n) --mhi;
		if(mhi < 0) mhi = -1;
		for( i= 0; i<=mhi; i++){
			s = 0.;
			for( h=0; h < m; h++)
				s += ((T)m_lpf_d[h])*x[2*i+h];
			y[i] = s;
		}

		/* fix up edge values */
		for( i= mhi+1; i<n2; i++){
			s = 0.;
			j = 2*i;
			for( h=0; h < m; h++){
				if(j >= n) j -= n;
				s += ((T)m_lpf_d[h])*x[j];
				j++;
			}
			y[i] = s;
		}
	}


	Matrix <T>
	idpwt2		(const Matrix <T> & wc, const int ell, const int J, T * temp)
	{

		Matrix <T> img (wc.DimVector());

		const int nr = wc.Height();
		const int nc = wc.Width();

		T *wcplo,*wcphi,*templo,*temphi,*temptop;
		int k,j,nj;
//		copydouble(wc,img,nr*nc);

		// assign dwt to result image
		img = wc;

		templo = &temp[nr];
		temphi = &temp[2*nr];
		temptop = &temp[3*nr];

		nj = 1;
		for( k=0; k < ell; k++) nj *=2;

		for( j=ell; j < J; j++){
			for( k=0; k < 2*nj; k++){
				unpackdouble(&img[0],nj,nc,k,templo);
				unpackdouble(&img[nj*nr],nj,nc,k,temphi);
				uplo(templo, nj,temp);
				uphi(temphi, nj, temptop);
				adddouble(temp,temptop,nj*2,temp);
				packdouble(temp,nj*2,nc,k,&img[0]);
			}

			for( k=0; k < 2*nj; k++){
				wcplo = &img[k*nr];
				wcphi = &img[k*nr + nj];
				copydouble(wcplo,temp,nj);
				uplo(wcplo, nj, templo);
				uphi(wcphi, nj, temphi);
				adddouble(templo,temphi,nj*2,wcplo);
			}
			nj *= 2;
		}

		return img;

	}


	void
	adddouble		(const T * const x, const T * const y, const int n, T * const z)
	{

		// TODO: use operator+ instead ...
		for (int i = n; i > 0; --i)
		{
			z[i] = x[i] + y[i];
		}

	}


	void
	uplo		(const T * const x, const int n, T * const y)
	{
	           int  meven, modd, i, h, j, mmax;
			   T s;

			   const int m = m_lpf_r.Dim (0);

	           /*lowpass version */

				/* away from edges */
		       meven = (m+1)/2; modd = m/2;
	           for( i= meven; i<n; i++){
			   		s = 0.;
					for( h=0; h < meven; h++)
						s += ((T)m_lpf_r[2*h])*x[i-h];
					y[2*i] = s;
					s = 0.;
					for( h=0; h < modd; h++)
						s += ((T)m_lpf_r[2*h+1])*x[i-h];
					y[2*i+1] = s;
				}

				/* fix up edge values */
				mmax = meven;
				if(mmax > n) mmax = n;
				for( i= 0; i < mmax; i++){

					s = 0.;
					j = i;
					for( h=0; h < meven; h++){
						if(j < 0) j += n;
						s += ((T)m_lpf_r[2*h])*x[j];
						--j;
					}
					y[2*i] = s;

					s = 0.;
					j = i;
					for( h=0; h < modd; h++){
						if(j < 0) j += n;
						s += ((T)m_lpf_r[2*h+1])*x[j];
						--j;
					}
					y[2*i+1] = s;
				}
	}


	void
	uphi		(const T * const x, const int n, T * const y)
	{
	           int  meven, modd, i, h, j, mmin;
			   T s;

			   const int m = m_hpf_r.Dim (0);

			   /*hipass version */
		       meven = (m+1)/2;
			   modd = m/2;

				/* away from edges */
	           for( i= 0; i+meven<n; i++){
			   		s = 0.;
					for( h=0; h < meven; h++)
						s += ((T)m_hpf_r[2*h])*x[i+h];
					y[2*i+1] = s;
					s = 0.;
					for( h=0; h < modd; h++)
						s += ((T)m_hpf_r[2*h+1])*x[i+h];
					y[2*i] = s;
				}

				/* fix up edge values */
				mmin = n-meven;
				if(mmin < 0) mmin = 0;
				for( i= mmin; i<n; i++){

					s = 0.;
					j = i;
					for( h=0; h < meven; h++){
						if(j >= n) j -= n;
						s += ((T)m_hpf_r[2*h])*x[j];
						j++;
					}
					y[2*i+1] = s;

					s = 0.;
					j = i;
					for( h=0; h < modd; h++){
						if(j >= n) j -= n;
						s += ((T)m_hpf_r[2*h+1])*x[j];
						j++;
					}
					y[2*i] = s;
				}
	}



};


# endif // __DWT2_HPP__
