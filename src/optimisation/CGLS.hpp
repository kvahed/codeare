/*
 * CGLS.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */

#ifndef SRC_OPTIMISATION_CGLS_HPP_
#define SRC_OPTIMISATION_CGLS_HPP_

#include <Linear.hpp>

namespace codeare {
namespace optimisation {

template<class T>
class CGLS : public Linear<T> {
public:
	CGLS (Optimisable* opt = 0, const int& verbosity = 0) :
		Linear<T>::Linear(opt, verbosity), _nrows(1), _ncols(1), _epsilon(1.e-15) {
	}
	virtual ~CGLS () {}
	int Solve () {

		int ret = 0;

		/*
        T rn, rno, xn, ts;
		Matrix<CT> p, r, x, q;
		vector<T> res;
        Vector<Matrix<cxfl> > vc;

        typedef typename Vector<CT>::iterator it_type;

		p = EH (m, sens, m_nx, m_fts)* m_ic;
		if (m_cgiter == 0)
			return p;

		r  = p;
        xn = real(p.dotc(p));
        rn = xn;

        if (m_verbose)
            vc.push_back (p);

		for (size_t i = 0; i < m_cgiter; i++) {
			res.push_back(rn/xn);
			if (i==0)
				x  = zeros<CT>(size(p));
			if (boost::math::isnan(res.at(i)) || res.at(i) <= m_cgeps)
				break;
			if (m_verbose)
				printf ("    %03lu %.7f\n", i, res.at(i));
			q  = EH(E(p * m_ic, sens, m_nx, m_fts), sens, m_nx, m_fts) * m_ic;

			if (m_lambda)
				q  += m_lambda * p;
			ts  = rn / real(p.dotc(q));
			x  += ts * p;
			r  -= ts * q;
			rno = rn;
			rn  = real(r.dotc(r));
			p  *= rn / rno;
			p  += r;
			if (m_verbose)
				vc.push_back(x * m_ic);
		}

        if (m_verbose) { // Keep intermediate results
            size_t cpsz = numel(x);
            x = Matrix<CT> (size(x,0), size(x,1), (m_nx[0] == 3) ? size(x,2) : 1, vc.size());
            it_type it = x.Begin();
            for (size_t i = 0; i < vc.size(); i++) {
                std::copy (vc[i].Begin(), vc[i].End(), it);
                it += cpsz;
            }
            vc.Clear();
        } else
            x *= m_ic;

 */
		return ret;
	}
protected:
	size_t _nrows, _ncols;
	T _epsilon, _rel_mat_err, _rel_rhs_err, _lambda;
};
}}

#endif /* SRC_OPTIMISATION_CGLS_HPP_ */
