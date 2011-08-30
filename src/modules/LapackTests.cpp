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

#include "LapackTests.hpp"


using namespace RRStrategy;


RRSModule::error_code
LapackTests::Process     () { 

	// SGEEV, DGEEV, CGEEV: Eigen value computation

	Matrix<double> hev;
	Matrix<double> hlev;
	Matrix<double> hrev;
	
	std::cout << "Testing EIG (dgeev) for helper:  ";
	std::cout << "INFO: " << m_real.begin()->second->EIG (&hev, &hlev, &hrev) << std::endl;

	Matrix<cplx> rev;
	Matrix<cplx> rlev;
	Matrix<cplx> rrev;

	std::cout << "Testing EIG (cgeev) for raw:     ";
	std::cout << "INFO: " << m_cplx.begin()->second->EIG (&rev, &rlev, &rrev)    << std::endl;

	// SEIG, DEIG, CEIG: Eigen value computation

	Matrix<double> dlsv;
	Matrix<double> drsv;
	Matrix<double> dsv;

	std::cout << "Testing SVD (dgesdd) for helper: ";
	std::cout << "INFO: " << m_real.begin()->second->SVD (&dlsv, &drsv, &dsv) << std::endl;

	Matrix<cplx>    rlsv;
	Matrix<cplx>    rrsv;
	Matrix<cplx>     rsv;

	std::cout << "Testing SVD (cgesdd) for raw:    ";
	std::cout << "INFO: " << m_cplx.begin()->second->SVD (&rlsv, &rrsv, &rsv)    << std::endl;

	Matrix<double> test;
	test.Dim(0) = 3;
	test.Dim(1) = m_real.begin()->second->Dim(0);
	test.Reset();
	test.Random();

	std::cout << "Testing mm multiplication (dgemm) for helper:     ";

	std::cout << test << std::endl;
	printf ("\n");
	test = test.prod((*m_real.begin()->second));
	std::cout << test << std::endl;

	m_cplx.begin()->second->Inv();

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new LapackTests;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
	delete p;
}
