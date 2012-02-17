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

#ifdef HAVE_INSIGHT
/*
#include "itkImage.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkOrientedImage.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
*/
template<class T> void
Matrix<T>::Resample (const Matrix<double>& f, const InterpMethod& im) {
	/*	
	typedef typename itk::OrientedImage< T, 3 > InputImageType;
	typedef typename itk::OrientedImage< T, 3 > OutputImageType;
	typedef typename itk::IdentityTransform< double, 3 > TransformType;
	typedef typename itk::LinearInterpolateImageFunction< InputImageType, double > InterpolatorType;
	typedef typename itk::ResampleImageFilter< InputImageType, InputImageType > ResampleFilterType;
	
	typename InterpolatorType::Pointer linterp = InterpolatorType::New();
	
	std::cout << "No interpolation attempted. Interpolation method unknown!" << std::endl;
	
	TransformType::Pointer trafo = TransformType::New();
	trafo->SetIdentity();
	
	typename InputImageType::SpacingType space;
	space[0] = 1.0/f[0];
	space[1] = 1.0/f[1];
	space[2] = 1.0/f[2];
	
	typedef typename InputImageType::SizeType::SizeValueType SizeValueType;
	typename InputImageType::SizeType size; 
	size[0] = static_cast<SizeValueType>(this->Dim(0));
	size[1] = static_cast<SizeValueType>(this->Dim(1));
	size[2] = static_cast<SizeValueType>(this->Dim(2));
	
	typename itk::OrientedImage< T, 3 >::Pointer input = itk::OrientedImage< T, 3 >::New();
	typename itk::OrientedImage< T, 3 >::Pointer output = itk::OrientedImage< T, 3 >::New();
	
	typename itk::Image< T, 3 >::IndexType ipos;
	ipos[0] = 0; ipos[1] = 0; ipos[2] = 0;
	typename itk::Image< T, 3 >::IndexType opos;
	opos[0] = 0; opos[1] = 0; opos[2] = 0;
	
	typename itk::Image< T, 3 >::RegionType ireg;
	ireg.SetSize(size);
	ireg.SetIndex(ipos);
	input->SetRegions(ireg);
	input->Allocate();
	
	typename itk::Image< T, 3 >::RegionType oreg;
	oreg.SetSize(size);
	ireg.SetIndex(opos);
	output->SetRegions(oreg);
	output->Allocate();
	
	for (size_t z = 0; z < this->Dim(2); z++)
		for (size_t y = 0; y < this->Dim(1); y++)
			for (size_t x = 0; x < this->Dim(0); x++) {
				ipos[0] = x; ipos[1] = y; ipos[2] = z;
				input->SetPixel (ipos, this->At(x,y,z));
			}
	
	typename ResampleFilterType::Pointer rs = ResampleFilterType::New();
	rs->SetInput( input );
	rs->SetTransform( trafo );
	rs->SetInterpolator( linterp );
	rs->SetOutputOrigin ( input->GetOrigin());
	rs->SetOutputSpacing ( space );
	rs->SetOutputDirection ( input->GetDirection());
	rs->SetSize ( size );
	rs->Update ();
	
	output = rs->GetOutput();
	
	this->Dim(0) = this->Dim(0)*f[0];
	this->Dim(1) = this->Dim(1)*f[1];
	this->Dim(2) = this->Dim(2)*f[2];
	this->Reset();

	this->Res(0) = this->Res(0)/f[0];
	this->Res(1) = this->Dim(1)/f[1];
	this->Res(2) = this->Dim(2)/f[2];
	
	for (size_t z = 0; z < this->Dim(2); z++)
		for (size_t y = 0; y < this->Dim(1); y++)
			for (size_t x = 0; x < this->Dim(0); x++) {
				opos[0] = x; opos[1] = y; opos[2] = z;
				this->At(x,y,z) = output->GetPixel (opos);
			}
	*/	
}



template<class T> Matrix<T>
Matrix<T>::Resample (const Matrix<double>& f, const InterpMethod& im) const {

	Matrix<T> res = (*this);
	return res->Resample(f, im);

}



template<class T> void
Matrix<T>::Resample (const double& f, const InterpMethod& im) {

	Matrix<double> mf (3,1);
	mf[0] = f; 
	mf[1] = f; 
	mf[2] = f; 

	return this->Resample(mf, im);

}


template<class T> Matrix<T>
Matrix<T>::Resample (const double& f, const InterpMethod& im) const {

	Matrix<T> res = (*this);

	Matrix<double> mf (3,1);
	mf[0] = f; 
	mf[1] = f; 
	mf[2] = f; 

	return res->Resample(mf, im);

}

#endif
