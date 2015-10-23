
/**
   Implementation of the loess algorithm in C++
   original code R project: http://www.r-project.org/
   http://svn.r-project.org/R/trunk/src/library/stats/src/lowess.c

   TEST:

   Using #RStats

   > T<-read.table('data.txt')
   > T
   V1   V2
   1  1  1.5
   2  3  2.0
   3  6  7.0
   4  9  8.0
   5 12 15.0

   > lowess(T$V1,T$V2)
   $x
   [1]  1  3  6  9 12

   $y
   [1] 1.5 2.0 7.0 8.0 8.0

   The same dataset, but using this class


   $ ./a.out < data.txt
   1	1.5	1.5
   3	2	2
   6	7	7
   9	8	8
   12	15	8

*/
#ifndef LOESS_STATS_H
#define LOESS_STATS_H
#include <stdint.h>
#include <Matrix.hpp>
#include <memory>

#define imin2(a,b) std::min(a,b)
#define imax2(a,b) std::max(a,b)
#define fmax2(a,b) std::max(a,b)

#define rPsort(x,n,k) std::partial_sort(&x[0],&x[n],&x[k]);

template<class T> class Loess {
	typedef typename TypeTraits<T>::RT RT;

public:
	Loess (const RT& span = 2./3., const size_t& steps = 1, const RT& speed = -1.,
			bool paranoid = false) : _span(span), _nsteps(steps), _speed(speed), _paranoid(paranoid) {};

	~Loess() {};

	RT _span; // span
	size_t _nsteps;       // robustness iterations
	RT _speed;   // speedup
	bool _paranoid;        // paranoid checks

	inline Matrix<T> lowess (const Matrix<T>& y) const {

		size_t n = y.Dim(0), howmany = numel(y)/n;
		assert(n>1);

		const Matrix<T>& x = linspace(0.,100.,n);

	    //TODO: Check for NaN
	    Matrix<T> ret(size(y));

	    RT delta=_speed;
	    if (delta < 0.0)
	    	delta = .01*(x[n-1]-x[0]);

#pragma omp parallel default (shared)
	    {
		    Vector<RT> rw(n);
		    Vector<RT> res(n);
#pragma omp	for schedule (static,1)
			for (size_t i = 0; i < howmany; ++i)
				clowess(&x[0], &y[0+i*n], n, _span, _nsteps, delta, &ret[0+i*n], &rw[0], &res[0]);
	    }
	    return ret;

    }


private:
	inline static T fcube(T x) { return x*x*x; }
	inline static T fsquare(T x) {return x*x; }
    
	void lowest (const RT *x, const RT *y, int n, const RT *xs, RT *ys, int nleft,
                 int nright, RT *w, bool userw, RT *rw, bool *ok) const {

        int nrt, j;
        RT a, b, c, h, h1, h9, r, range;

        x--;
        y--;
        w--;
        rw--;

        range = x[n]-x[1];
        h = fmax2(*xs-x[nleft], x[nright]-*xs);
        h9 = 0.999*h;
        h1 = 0.001*h;

        /* sum of weights */

        a = 0.;
        j = nleft;
        while (j <= n) {

			/* compute weights */
			/* (pick up all ties on right) */

			w[j] = 0.;
			r = fabs(x[j] - *xs);
			if (r <= h9) {
			    if (r <= h1)
                    w[j] = 1.;
			    else
                    w[j] = fcube(1.-fcube(r/h));
			    if (userw)
                    w[j] *= rw[j];
			    a += w[j];
			}
			else if (x[j] > *xs)
			    break;
			j = j+1;
        }

        /* rightmost pt (may be greater */
        /* than nright because of ties) */

        nrt = j-1;
        if (a <= 0.)
			*ok = false;
        else {
			*ok = true;

			/* weighted least squares */
			/* make sum of w[j] == 1 */

			for(j=nleft ; j<=nrt ; j++)
			    w[j] /= a;
			if (h > 0.) {
			    a = 0.;

			    /*  use linear fit */
			    /* weighted center of x values */

			    for(j=nleft ; j<=nrt ; j++)
                    a += w[j] * x[j];
			    b = *xs - a;
			    c = 0.;
			    for(j=nleft ; j<=nrt ; j++)
                    c += w[j]*fsquare(x[j]-a);
			    if (sqrt(c) > 0.001*range) {
                    b /= c;

                    /* points are spread out */
                    /* enough to compute slope */

                    for(j=nleft; j <= nrt; j++)
                        w[j] *= (b*(x[j]-a) + 1.);
			    }
			}
			*ys = 0.;
			for(j=nleft; j <= nrt; j++)
			    *ys += w[j] * y[j];
        }
    }

	void clowess (const RT  *x, const RT *y, int n, RT f, int nsteps, RT delta, RT *ys, RT *rw,
			RT *res) const {

        int i, iter, j, last, m1, m2, nleft, nright, ns;
        bool ok;
        RT alpha, c1, c9, cmad, cut, d1, d2, denom, r, sc;

        if (n < 2) {
            ys[0] = y[0]; return;
        }

        /* nleft, nright, last, etc. must all be shifted to get rid of these: */
        x--;
        y--;
        ys--;

        /* at least two, at most n points */
        ns = imax2(2, imin2(n, (int)(f*n + 1e-7)));
#ifdef DEBUG_lowess
        REprintf("lowess(): ns = %d\n", ns);
#endif

        /* robustness iterations */

        iter = 1;
        while (iter <= nsteps+1) {
            nleft = 1;
            nright = ns;
            last = 0;	/* index of prev estimated point */
            i = 1;		/* index of current point */

            for(;;) {
                if (nright < n) {

                    /* move nleft,  nright to right */
                    /* if radius decreases */

                    d1 = x[i] - x[nleft];
                    d2 = x[nright+1] - x[i];

                    /* if d1 <= d2 with */
                    /* x[nright+1] == x[nright], */
                    /* lowest fixes */

                    if (d1 > d2) {

                        /* radius will not */
                        /* decrease by */
                        /* move right */

                        nleft++;
                        nright++;
                        continue;
                    }
                }

                /* fitted value at x[i] */

                lowest(&x[1], &y[1], n, &x[i], &ys[i],
                       nleft, nright, res, iter>1, rw, &ok);
                if (!ok) ys[i] = y[i];

                /* all weights zero */
                /* copy over value (all rw==0) */

                if (last < i-1) {
                    denom = x[i]-x[last];

                    /* skipped points -- interpolate */
                    /* non-zero - proof? */

                    for(j = last+1; j < i; j++) {
                        alpha = (x[j]-x[last])/denom;
                        ys[j] = alpha*ys[i] + (1.-alpha)*ys[last];
                    }
                }

                /* last point actually estimated */
                last = i;

                /* x coord of close points */
                cut = x[last]+delta;
                for (i = last+1; i <= n; i++) {
                    if (x[i] > cut)
                        break;
                    if (x[i] == x[last]) {
                        ys[i] = ys[last];
                        last = i;
                    }
                }
                i = imax2(last+1, i-1);
                if (last >= n)
                    break;
            }
            /* residuals */
            for(i = 0; i < n; i++)
                res[i] = y[i+1] - ys[i+1];

            /* overall scale estimate */
            sc = 0.;
            for(i = 0; i < n; i++) sc += fabs(res[i]);
            sc /= n;

            /* compute robustness weights */
            /* except last time */

            if (iter > nsteps)
                break;
            /* Note: The following code, biweight_{6 MAD|Ri|}
               is also used in stl(), loess and several other places.
               --> should provide API here (MM) */
            for(i = 0 ; i < n ; i++)
                rw[i] = fabs(res[i]);

            /* Compute   cmad := 6 * median(rw[], n)  ---- */
            /* FIXME: We need C API in R for Median ! */
            m1 = n/2;
            /* partial sort, for m1 & m2 */
            rPsort(rw, n, m1);
            if(n % 2 == 0) {
                m2 = n-m1-1;
                rPsort(rw, n, m2);
                cmad = 3.*(rw[m1]+rw[m2]);
            }
            else { /* n odd */
                cmad = 6.*rw[m1];
            }

            if(cmad < 1e-7 * sc) /* effectively zero */
                break;
            c9 = 0.999*cmad;
            c1 = 0.001*cmad;
            for(i = 0 ; i < n ; i++) {
                r = fabs(res[i]);
                if (r <= c1)
                    rw[i] = 1.;
                else if (r <= c9)
                    rw[i] = fsquare(1.-fsquare(r/cmad));
                else
                    rw[i] = 0.;
            }
            iter++;
        }
    }
};

#endif

