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

#include "Lapack.hpp"

class MCGLS {

public:
	
	template<class T> static Matrix<T> 
	Pinv (Matrix<T>& A, Matrix<T>& b, const size_t& maxit, const double& conv, const double& lambda) {

		size_t ah = A.Height();
		size_t aw = A.Width();
		size_t bh = b.Height();
		size_t bw = b.Width();

		assert (bw == 1);  // Column vector x
		assert (ah == bh); // Check inner dimensions of A'*x. 

		ticks tic = getticks();

		Matrix<T> p = Lapack::GEMM (A, b, 'C');
		Matrix<T> r = p;

		Matrix<T> x (p.Dim()); 
		Matrix<T> q;

		T         ts;

		float     rn = 0.0;
		float     xn = pow(p.Norm().real(), 2);

		std::vector<double> res;

		for (size_t i = 0; i < maxit; i++) {
			
			rn  = pow(r.Norm().real(), 2);
			res.push_back(rn/xn);

			if (std::isnan(res.at(i)) || res.at(i) <= conv) break;

			if (i % 5 == 0 && i > 0)                        printf ("\n");
			printf ("    %03lu %.7f", i, res.at(i));
			
			q  = Lapack::GEMM (A, p);
			q  = Lapack::GEMM (A, q, 'C');
			q += (T(lambda)*p);
			
			ts  = (rn / (p.dotc(q)));
			x  += (p * ts);
			r  -= (q * ts);
			p  *= cxfl(pow(r.Norm().real(), 2)/rn);
			p  += r;
			
		}
		
		printf ("\n  MCGLS time: %.4f s\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
		
		return x;
		
	}

};
