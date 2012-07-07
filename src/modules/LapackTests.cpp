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

#include "LapackTests.hpp"
#include "Lapack.hpp"
#include "MCGLS.hpp"
#include "CX.hpp"
#include "IO.hpp"
#include "Statistics.hpp"
#include "Creators.hpp"

using namespace RRStrategy;


RRSModule::error_code
LapackTests::Process     () { 

	// SGEEV, DGEEV, CGEEV, ZGEEV: Eigen value computation

	std::cout << "General LA tests ----------- \n";

	Matrix<cxfl>&   cf = GetCXFL ("cf");
	Matrix<cxdb>    cd = (Matrix<cxdb>) cf;
	Matrix<double>& rd = GetRLDB ("rd");

	std::cout << "Testing eig (cgeev) ----------- \n";

	Matrix<cxdb>   hev (rd.Dim(0), 1);
	Matrix<double> hlev(rd.Dim(0));
	Matrix<double> hrev(rd.Dim(0));
	Matrix<double> dtmp = rd;
	
	std::cout << "RLDB IN: \n";
	std::cout << rd << std::endl;

	eig (dtmp, hev, hlev, hrev, 'V', 'V');
	
	std::cout << "OUT: \n" <<  hev;
 	std::cout << "------------------------------- \n\n";


	std::cout << "Testing eig (cgeev) ----------- \n";

	Matrix<cxfl> rev  (cf.Dim(0), 1);
	Matrix<cxfl> rlev (cf.Dim(0));
	Matrix<cxfl> rrev (cf.Dim(0));
	Matrix<cxfl> rtmp = rd;

	std::cout << "CXFL IN: \n";
	std::cout << cf << std::endl;

	eig (rtmp, rev, rlev, rrev);

	std::cout << "OUT: \n" <<  rev;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing svd (dgesdd) ---------- \n";

	Matrix<double> Ar = rand<double> (8,3);
	Matrix<double> du;
	Matrix<double> dv;
	Matrix<double> ds;

	std::cout << "RLDB IN: \n";
	std::cout << Ar << std::endl;
	svd (Ar, ds, du, dv, 'O');
	std::cout << "OUT s:\n" <<  ds  << std::endl;
	std::cout << "OUT u:\n" <<  du  << std::endl;
	std::cout << "OUT v:\n" <<  dv;
	std::cout << "------------------------------- \n\n";
	
	std::cout << "Testing svd (dgesdd) ---------- \n";

	Matrix<cxfl> Ac = rand<cxfl> (3,8); 

	Matrix<cxfl>  cu;
	Matrix<cxfl>  cv;
	Matrix<float> cs;

	std::cout << "CXFL IN: \n";
	std::cout << Ac << std::endl;
	svd (Ac, cs, cu, cv, 'S');
	std::cout << "OUT s:\n" <<  cs  << std::endl;
	std::cout << "OUT u:\n" <<  cu  << std::endl;
	std::cout << "OUT v:\n" <<  cv;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing Inv (dgetri / dgetrf) - \n";

	dtmp = rand<double> (5,5);
	std::cout << "RLDB IN: \n";
	std::cout << dtmp << std::endl;
	dtmp = inv(dtmp);
	std::cout << "OUT:\n" <<  dtmp;
	std::cout << "------------------------------- \n\n";
	
	std::cout << "Testing Inv (dgetri / dgetrf) - \n";

	Matrix<cxdb>ctmp = rand<cxdb> (5,5);
	std::cout << "CXFL IN: \n";
	std::cout << ctmp << std::endl;
	ctmp = inv(ctmp);
	std::cout << "OUT:\n" <<  ctmp;
	std::cout << "------------------------------- \n\n";
	
	std::cout << "Testing Pinv (dgelsd) ---------- \n";

	Matrix<double> da = rand<double> (5,3);

	std::cout << "RLDB IN: \n";
	std::cout << da << std::endl;
	da = pinv (da);

	std::cout << "OUT inv:\n" <<  da;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing Pinv (zgelsd) ---------- \n";

	Matrix<cxfl> db = rand<cxfl> (3,5);

	std::cout << "CXDB IN: \n";
	std::cout << db << std::endl;
	db = pinv (db);

	std::cout << "OUT inv:\n" <<  db;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing mm multipl. (dgemm) --- \n";

	Matrix<cxfl> dc = rand<cxfl> (3,8);
	
	std::cout << "A = [\n";
	std::cout << db;
	std::cout << "];" << std::endl;
	std::cout << "B = [\n";
	std::cout << dc;
	std::cout << "];" << std::endl;

	printf ("\n");
	Matrix<cxfl> dd = gemm (db, dc);
	std::cout << "A*B:\n" <<  dd;
	dd = gemm (dc, db, 'T', 'T');
	std::cout << "B.'*A.' :\n" <<  dd;
	dd = gemm (dc, db, 'C', 'C');
	std::cout << "B'*A' :\n" <<  dd;
	db = !db;

	std::cout << "A  :\n" <<  db;
	std::cout << "B' :\n" <<  dc;

	dd = gemm (db, dc, 'C');
	std::cout << "A*(B').' :\n" <<  dd;
	std::cout << "------------------------------- \n\n";

	Matrix<cxfl> A = rand<cxfl> (3,8); 

	Matrix<cxfl> b = rand<cxfl> (A.Height());
	
	std::cout << b.dotc(b) << std::endl;

	std::cout << "A  :\n" <<  A;
	std::cout << "x' :\n" <<  b;

	ticks tic = getticks();

	tic = getticks();

	Matrix<cxfl> y;

	y = gemm (A, b, 'C');
	printf ("gemm. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	std::cout << "y' :\n" <<  y;
	
	y  = gemv (A, b, 'C');
	printf ("gemv. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	std::cout << "y' :\n" <<  y;


	A = rand<cxfl> (1000,8);

	A = cov (A);
	std::cout << "cov (A)  :\n" <<  A;
	A = chol (A);
	std::cout << "chol(cov(A)) :\n" <<  A;

	A = GetCXFL("EM");
	b = GetCXFL("PA");

	Matrix<cxfl> x = AddMatrix ("x", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>  (1)));
	x = MCGLS (A, b, 100, 5.0e-6, 1.0e-3);

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new LapackTests;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
	delete p;
}
