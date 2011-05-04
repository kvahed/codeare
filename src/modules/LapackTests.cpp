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

	Matrix<raw> ev;
	Matrix<double> hlev;
	Matrix<double> hrev;
	
	std::cout << "Testing EIG (dgeev) for helper:  ";
	std::cout << "INFO: " << m_helper.EIG (true, &ev, &hlev, &hrev) << std::endl;

	Matrix<raw> rlev;
	Matrix<raw> rrev;

	std::cout << "Testing EIG (cgeev) for raw:     ";
	std::cout << "INFO: " << m_raw.EIG (true, &ev, &rlev, &rrev)    << std::endl;

	// SEIG, DEIG, CEIG: Eigen value computation

	Matrix<double> dlsv;
	Matrix<double> drsv;
	Matrix<double> dsv;

	std::cout << "Testing SVD (dgesdd) for helper: ";
	std::cout << "INFO: " << m_helper.SVD (true, &dlsv, &drsv, &dsv) << std::endl;

	Matrix<raw>    rlsv;
	Matrix<raw>    rrsv;

	std::cout << "Testing SVD (cgesdd) for raw:    ";
	std::cout << "INFO: " << m_raw.SVD (true, &rlsv, &rrsv, &dsv)    << std::endl;

	Matrix<short> test;
	test.Dim(COL) = m_pixel.Dim(LIN);
	test.Dim(LIN) = 3;
	test.Reset();
	test.Random();

	std::cout << "Testing mm multiplication (dgemm) for helper:     ";
	Matrix<short> product = test.prod(m_pixel);
	
	std::cout << test;
	printf ("\n");
	std::cout << product;

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new LapackTests;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
	delete p;
}
