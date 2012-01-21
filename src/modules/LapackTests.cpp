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

	Matrix<cxfl>& cf = GetCXFL ("cf");
	Matrix<double>& rd = GetRLDB ("rd");

	Matrix<double> hev;
	Matrix<double> hlev;
	Matrix<double> hrev;
	
	std::cout << "Testing EIG (dgeev) for helper:  ";
	std::cout << "INFO: " << rd.EIG (&hev, &hlev, &hrev) << std::endl;

	Matrix<cxfl> rev;
	Matrix<cxfl> rlev;
	Matrix<cxfl> rrev;

	std::cout << "Testing EIG (cgeev) for raw:     ";
	std::cout << "INFO: " << cf.EIG (&rev, &rlev, &rrev)    << std::endl;

	// SEIG, DEIG, CEIG: Eigen value computation

	Matrix<double> dlsv;
	Matrix<double> drsv;
	Matrix<double> dsv;

	std::cout << "Testing SVD (dgesdd) for helper: ";
	std::cout << "INFO: " << rd.SVD (&dlsv, &drsv, &dsv) << std::endl;

	Matrix<cxfl>    rlsv;
	Matrix<cxfl>    rrsv;
	Matrix<cxfl>     rsv;

	std::cout << "Testing SVD (cgesdd) for raw:    ";
	std::cout << "INFO: " << cf.SVD (&rlsv, &rrsv, &rsv)    << std::endl;

	Matrix<double> test;
	test.Dim(0) = 3;
	test.Dim(1) = rd.Dim(0);
	test.Reset();
	test.Random();

	std::cout << "Testing mm multiplication (dgemm) for helper:     ";

	std::cout << test << std::endl;
	printf ("\n");
	test = test.prod(rd);
	std::cout << test << std::endl;

	cf.Inv();

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new LapackTests;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
	delete p;
}
