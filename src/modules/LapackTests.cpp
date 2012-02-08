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
#include "Lapack.hpp"
#include "MCGLS.hpp"
#include "CX.hpp"

using namespace RRStrategy;


RRSModule::error_code
LapackTests::Process     () { 

	// SGEEV, DGEEV, CGEEV, ZGEEV: Eigen value computation


	Matrix<cxfl>&   cf = GetCXFL ("cf");
	Matrix<double>& rd = GetRLDB ("rd");

	std::cout << "Testing EIG (cgeev) ----------- \n";

	Matrix<cxdb>   hev (rd.Dim(0), 1);
	Matrix<double> hlev(rd.Dim(0));
	Matrix<double> hrev(rd.Dim(0));
	Matrix<double> dtmp = rd;
	
	std::cout << "RLDB IN: \n";
	std::cout << rd << std::endl;

	Lapack::EIG (dtmp, hev, hlev, hrev, 'V', 'V');
	
	std::cout << "OUT: \n" <<  hev;
	std::cout << "------------------------------- \n\n";


	std::cout << "Testing EIG (cgeev) ----------- \n";

	Matrix<cxfl> rev  (cf.Dim(0), 1);
	Matrix<cxfl> rlev (cf.Dim(0));
	Matrix<cxfl> rrev (cf.Dim(0));
	Matrix<cxfl> rtmp = rd;

	std::cout << "CXFL IN: \n";
	std::cout << cf << std::endl;

	Lapack::EIG (rtmp, rev, rlev, rrev);

	std::cout << "OUT: \n" <<  rev;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing SVD (dgesdd) ---------- \n";

	Matrix<double> du (MIN(rd.Height(),rd.Width()));
	Matrix<double> dv (MAX(rd.Height(),rd.Width()));
	Matrix<double> ds (MIN(rd.Height(),rd.Width()), 1);

	std::cout << "RLDB IN: \n";
	std::cout << rd << std::endl;
	std::cout << "INFO: " << Lapack::SVD (rd, ds, du, dv, 'A') << "\n" << std::endl;
	std::cout << "OUT s:\n" <<  ds  << std::endl;
	std::cout << "OUT u:\n" <<  du  << std::endl;
	std::cout << "OUT v:\n" <<  dv;
	std::cout << "------------------------------- \n\n";
	
	std::cout << "Testing SVD (dgesdd) ---------- \n";

	Matrix<cxfl>  cu (MIN(cf.Height(),cf.Width()));
	Matrix<cxfl>  cv (MAX(cf.Height(),cf.Width()));
	Matrix<float> cs (MIN(cf.Height(),cf.Width()), 1);

	std::cout << "RLDB IN: \n";
	std::cout << cf << std::endl;
	std::cout << "INFO: " << Lapack::SVD (cf, cs, cu, cv, 'A') << "\n" << std::endl;
	std::cout << "OUT s:\n" <<  cs  << std::endl;
	std::cout << "OUT u:\n" <<  cu  << std::endl;
	std::cout << "OUT v:\n" <<  cv;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing Inv (dgetri / dgetrf) - \n";

	dtmp = Matrix<double>::Random2D (5);
	std::cout << "RLDB IN: \n";
	std::cout << dtmp << std::endl;
	dtmp = Lapack::Inv(dtmp);
	std::cout << "OUT:\n" <<  dtmp;
	std::cout << "------------------------------- \n\n";
	
	std::cout << "Testing Inv (dgetri / dgetrf) - \n";

	Matrix<cxdb>ctmp = Matrix<cxdb>::Random2D (5);
	std::cout << "CXFL IN: \n";
	std::cout << ctmp << std::endl;
	ctmp = Lapack::Inv(ctmp);
	std::cout << "OUT:\n" <<  ctmp;
	std::cout << "------------------------------- \n\n";
	

	std::cout << "Testing Pinv (dgelsd) ---------- \n";

	Matrix<double> da (5,3);
	da.Random();

	std::cout << "RLDB IN: \n";
	std::cout << da << std::endl;
	da = Lapack::Pinv (da);

	std::cout << "OUT inv:\n" <<  da;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing Pinv (zgelsd) ---------- \n";

	Matrix<cxfl> db (3,5);
	db.Random();

	std::cout << "CXDB IN: \n";
	std::cout << db << std::endl;
	db = Lapack::Pinv (db);

	std::cout << "OUT inv:\n" <<  db;
	std::cout << "------------------------------- \n\n";

	std::cout << "Testing mm multipl. (dgemm) --- \n";

	Matrix<cxfl> dc (3,4) ;
	dc.Random();
	
	std::cout << "A = [\n";
	std::cout << db;
	std::cout << "];" << std::endl;
	std::cout << "B = [\n";
	std::cout << dc;
	std::cout << "];" << std::endl;

	printf ("\n");
	Matrix<cxfl> dd = Lapack::GEMM (db, dc);
	std::cout << "A*B:\n" <<  dd;
	dd = Lapack::GEMM (dc, db, 'T', 'T');
	std::cout << "B.'*A.' :\n" <<  dd;
	dd = Lapack::GEMM (dc, db, 'C', 'C');
	std::cout << "B'*A' :\n" <<  dd;
	db = !db;

	std::cout << "A  :\n" <<  db;
	std::cout << "B' :\n" <<  dc;

	dd = Lapack::GEMM (db, dc, 'C', 'N');
	std::cout << "A*(B').' :\n" <<  dd;
	std::cout << "------------------------------- \n\n";

	Matrix<cxfl> x (dc.Width(),1);
	x.Random();

	dc = !dc;
	std::cout << "A  :\n" <<  dc;
	std::cout << "x' :\n" <<  x;

	Matrix<cxfl> y = Lapack::GEMV (dc, x, 'C');

	std::cout << "y' :\n" <<  y;



	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new LapackTests;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
	delete p;
}
