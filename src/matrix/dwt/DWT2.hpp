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
//    	m_lpf_d [0] = -0.0106; m_lpf_d [1] = 0.0329; m_lpf_d [2] = 0.0308; m_lpf_d [3] = -0.1870;
//    	m_lpf_d [4] = -0.0280; m_lpf_d [5] = 0.6309; m_lpf_d [6] = 0.7148; m_lpf_d [7] = 0.2304;
        m_lpf_r [0] = 0.2304; m_lpf_r [1] = 0.7148; m_lpf_r [2] = 0.6309; m_lpf_r [3] = -0.0280;
        m_lpf_r [4] = -0.1870; m_lpf_r [5] = 0.0308; m_lpf_r [6] = 0.0329; m_lpf_r [7] = -0.0106;

        m_lpf_d = m_lpf_r;

    	// set low pass filter (haar)
//    	m_hpf_d [0] = -1; m_hpf_d [1] = 1;
//    	m_hpf_d *= norm_factor;
//        m_hpf_r [0] = 1; m_hpf_r [1] = -1;
//        m_hpf_r *= norm_factor;
//    	m_hpf_d [0] = -0.2304; m_hpf_d [1] = 0.7148; m_hpf_d [2] = -0.6409; m_hpf_d [3] = -0.0280;
//    	m_hpf_d [4] = 0.1870; m_hpf_d [5] = 0.0308; m_hpf_d [6] = -0.0329; m_hpf_d [7] = -0.0106;
//        m_hpf_d [0] = -0.0106; m_hpf_d [1] = -0.0329; m_hpf_d [2] = 0.0308; m_hpf_d [3] = 0.1870;
//        m_hpf_d [4] = -0.0280; m_hpf_d [5] = -0.6309; m_hpf_d [6] = 0.7148; m_hpf_d [7] = -0.2304;
        mirrorfilt((double*)m_lpf_d.Memory(0), (double*)m_hpf_d.Memory(0),m_lpf_d.Dim(0));
        mirrorfilt((double*)m_lpf_r.Memory(0), (double*)m_hpf_r.Memory(0),m_lpf_r.Dim(0));

//    	m_hpf_r = m_hpf_d;

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
        	const int ell = 1;

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
        	const int ell = 1;

        	// call dpwt2
        	res = idpwt2 (m, ell, J, temp);

        	free (temp);

        }

		return res;

	}


	void
	mirrorfilt              (double * lpf, double * hpf, int length)
	{
	    int i,isign;
	    isign = 1;
	    for(i=0; i < length; i++){
	        *hpf++ = isign * *lpf++;
	        isign *= -1;
	    }
	}



	/**
	 *  COPIED FROM WAVELAB IMPLEMENTATION
	 */

	Matrix <T>
	dpwt2		(const Matrix <T> & sig, const int ell, const int J, T * temp)
	{

		Matrix <T> res (sig.DimVector ());

		T *wcplo,*wcphi,*templo,*temphi;
//		copydouble(&sig[0],&res[0],sig.Height()*sig.Width());

		// assign signal to result matrix
		res = sig;

		templo = &temp[sig.Height ()];
		temphi = &temp[2*sig.Height()];

		int num_rows = sig.Height();
		int num_cols = sig.Width();
		int side_length = sig.Height();

		// loop over levels of DWT
		for (int j = (J-1); j >= ell; --j)
		{

		    // loop over columns of image
			for (int col=0; col < side_length; col++)
			{

			    // access to lowpass part of DWT
				wcplo = &res[col*num_rows];
				// access to highpass part of DWT
				wcphi = &res[col*num_rows + side_length/2];

				// copy part of image to temp memory
				copydouble (wcplo, temp, side_length);

				// apply low pass filter on column and write to result matrix
				downlo (temp, side_length, wcplo);
				// apply high pass filter on column and wirte to result matrix
				downhi (temp, side_length, wcphi);

			} // loop over columns

			// loop over rows of image
			for (int row=0; row < side_length; row++)
			{

			    // copy row-th row of imag to temp
				unpackdouble (res.Memory(0), side_length, num_cols, row, temp);

				// apply low pass filter on row and write to temp mem
				downlo (temp, side_length, templo);
				// apply high pass filter on row and write to temp mem
				downhi (temp, side_length, temphi);

				// write temp lowpass result to result matrix
				packdouble (templo, side_length/2, num_cols, row, &res[0]);
				// write temp highpass result to result matrix
				packdouble (temphi, side_length/2, num_cols, row, &res[side_length/2*num_rows]);

			} // loop over rows of image

			// reduce dimension for next level
			side_length = side_length/2;

		} // loop over levels of DWT

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
	copydouble		(const T * const src, T * const dest, const int n)
	{

		// TODO: use operator= instead ...
		for (int i = 0; i < n; ++i)
		{
			dest [i] = src [i];
		}

	}


    /**
     * @brief               Perform highpass part of one level forward 1D DWT on signal with given side_length.
     *
     * @param  signal       Signal to be transformed.
     * @param  side_length  Side length of current DWT level.
     * @param  dwt_hgih     Resulting DWT.
     */
	void
	downhi			(const T * const signal, const int side_length, T * const dwt_high)
	{
		int n2, mlo, j;
		T s;

		int filter_length = m_hpf_d.Dim (0);

		/* highpass version */

		// half of side_length
		n2 = side_length/2;

		// half of filter length
		mlo = filter_length/2 - 1;

		// adjust lower bound if to low
		if (2*mlo+1 - (filter_length-1) < 0)
		    mlo++;

		// loop over pixels of dwt_high
		for (int i = mlo; i < n2; i++)
		{

		    // result of convolution
			s = 0.;

			// perform convolution
			for (int h = 0; h < filter_length; h++)
				s += ((T) m_hpf_d [h]) * signal [2*i+1-h];

			// assign result of convolution
			dwt_high[i] = s;

		} // loop over pixels of dwt_high

		// case: filter_length > side_length => only edge values
		if (mlo > n2)
		    mlo = n2;


		/* fix up edge values */

		// loop over edge values
		for (int i = 0; i < mlo; i++)
		{

		    // result of convolution
			s = 0.;

			// start signal index for convolution
			j = 2*i+1;

			// loop over filter elements
			for (int h = 0; h < filter_length; h++)
			{

			    // adjust index if it exceeds side_length
				if (j < 0)
				    j += side_length;

				s += ((T) m_hpf_d [h]) * signal [j];

				// update index
				--j;

			}

			// assign result of convolution
			dwt_high[i] = s;

		} // loop over edge values

	}


	/**
	 * @brief               Perform lowpass part of one level forward 1D DWT on signal with given side_length.
	 *
	 * @param  signal       Signal to be transformed.
	 * @param  side_length  Side length of current DWT level.
	 * @param  dwt_low      Resulting DWT.
	 */
	void
	downlo  		        (const T * const signal, const int side_length, T * const dwt_low)
	{

		int n2, mlo, mhi, j;
		T s;

		int filter_length = m_lpf_d.Dim (0);

		/*lowpass version */

		// half of side_length (length of dwt_low)
		n2 = side_length/2;

		// half of filter_length
		mlo = filter_length /2;

		// upper bound for "normal" convolution
		mhi = n2 - mlo;

		// upper bound to high
		if (2*mhi + (filter_length-1) >= side_length)
		    --mhi;
		// upper bound to low
		if (mhi < 0)
		    mhi = -1;

		// loop over pixels of dwt_low
		for (int i= 0; i<=mhi; i++)
		{

		    // result of convolution
			s = 0.;
			// apply low pass filter (convolution)
			for (int h = 0; h < filter_length; h++)
				s += ((T) m_lpf_d [h]) * signal [2*i+h];
			// assign result of convolution
			dwt_low [i] = s;

		} // loop over pixels of dwt_low


		/* fix up edge values */

		// loop over edge values (periodic boundary)
		for (int i = mhi+1; i < n2; i++)
		{

		    // result of convolution
			s = 0.;

			// start signal index for convolution
			j = 2*i;

			// loop over filter elements
			for (int h = 0; h < filter_length; h++){

			    // adjust index if it exceeds current side_length
				if (j >= side_length)
				    j -= side_length;
				s += ((T) m_lpf_d [h]) * signal [j];

				// update index
				j++;

			} // perform convolution

			// assign result of convolution
			dwt_low [i] = s;

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
		for (int i = 0; i < n; ++i)
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
