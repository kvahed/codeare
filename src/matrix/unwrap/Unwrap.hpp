/*
  3D Unwrapping algorithm
  Rhodri Cusack 2000-2006
  Algorithm described in 
  Cusack, R. & Papadakis, N. (2002) 
  "New robust 3-D phase unwrapping algorithms: application to magnetic field mapping and undistorting echoplanar images."
  Neuroimage. 2002 Jul;16(3 Pt 1):754-64.
  Distributed under GNU public license http://www.gnu.org/copyleft/gpl.html
  Comments to rhodri.cusack@mrc-cbu.cam.ac.uk
  Version 3.00 Oct 2006 adapted for matlab
*/

#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include <math.h>


size_t m_bsx,m_bsy,m_bsz;
container<unsigned short> flag;

struct FIELD {
    double d;	/* value of field */
    size_t p[3];  /* at this offset */
    int x,y,z;
};

struct QUEUEENTRY {
    int x,y,z;
    size_t p;
    double v;
};

struct FIELD_2 {
    double d;	/*value of field*/
    size_t p;		/*at this offset*/
    int x,y,z;
};



const static size_t NUMQUEUES = 10000;
const static size_t DEFAULT_NB = NUMQUEUES;
const static size_t DEFAULTBLOCK = 100;
const static size_t BLOCKINCREMENT = 500;

size_t m_bot[NUMQUEUES], m_top[NUMQUEUES];
size_t m_size[NUMQUEUES], m_chunk[NUMQUEUES], m_sizeb[NUMQUEUES];
char *m_q[NUMQUEUES];

inline static void
InitQueue (const int queuenum, const int chunk) {
    
    m_chunk[queuenum] = chunk;
    m_bot[queuenum]   = 0;
    m_top[queuenum]   = 0;
    m_size[queuenum]  = DEFAULTBLOCK;
    m_sizeb[queuenum] = m_size[queuenum]*m_chunk[queuenum];
    
    m_q[queuenum]=0;

}

inline static void
TerminateQueue(const int queuenum) {
    
    if (m_q[queuenum])
    	delete m_q[queuenum];

    m_q[queuenum]=0;

}


inline static int
Push(const int queuenum, const void *x) {
    
    if (!m_q[queuenum])
    	m_q[queuenum] = new char[m_sizeb[queuenum]];

    if (!m_q[queuenum]) {
        printf("Out of memory - could not generate new point queue.\n");
        return(-1);
    }
    
    int stacksize = m_top[queuenum]-m_bot[queuenum];
    if (stacksize<0)
    	stacksize+=m_size[queuenum];

    memcpy(m_q[queuenum]+m_top[queuenum]*m_chunk[queuenum],x,m_chunk[queuenum]);
    m_top[queuenum]++;
    
    if (m_top[queuenum]==m_size[queuenum])
    	m_top[queuenum]=0;
    
    if (m_top[queuenum]==m_bot[queuenum]) {
        
        char *newq;
        size_t newsize, newsizeb, abovebot_b;

        newsize=m_size[queuenum]+BLOCKINCREMENT;
        newsizeb=newsize*m_chunk[queuenum];
        newq=new char[newsizeb];
        if (!newq) {
            printf("Out of memory - point queue full.\n");
            return(-1);
        }
        
        abovebot_b=(m_size[queuenum]-m_bot[queuenum])*m_chunk[queuenum];
        if (abovebot_b>0)
        	memcpy(newq,(char *)m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],abovebot_b);
        // then, do m_top downwards
        if (m_top[queuenum]!=0)
        	memcpy((char *) newq+abovebot_b,m_q[queuenum],m_top[queuenum]*m_chunk[queuenum]);
        m_bot[queuenum] =0;
        m_top[queuenum] = m_size[queuenum];
        m_size[queuenum] = newsize;
        m_sizeb[queuenum] = newsizeb;
        delete m_q[queuenum]; // recover old memory
        m_q[queuenum]=newq;

    }
    
    return 0;
    
}

static inline int
Pop (const int queuenum, void *x) {
    
    if (m_bot[queuenum]==m_top[queuenum])
    	return(1);
    
    memcpy(x,m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],m_chunk[queuenum]);
    m_bot[queuenum]++;

    if (m_bot[queuenum]==m_size[queuenum])
    	m_bot[queuenum]=0;
    
    return 0;
}


template<class T> static inline void
Check(int primaryqueuenum, QUEUEENTRY* qe, size_t offp, int offx, int offy, int offz,
	  const container<size_t>& dim, const Matrix<T>& phase, Matrix<T>& unwrapped) {

    // first check bounds
    QUEUEENTRY nqe;

    nqe.x=qe->x+offx;
    if ((nqe.x<0)||(nqe.x>=dim[0]))
    	return;

    nqe.y=qe->y+offy;
    if ((nqe.y<0)||(nqe.y>=dim[1]))
    	return;

    nqe.z=qe->z+offz;
    if ((nqe.z<0)||(nqe.z>=dim[2]))
    	return;

    nqe.p=qe->p+offp;

    if (flag[nqe.p])
    	return; // Already been here

    // Actually do unwrap
    int wholepis = int((phase[nqe.p]-qe->v)/PI);
    
    if (wholepis>=1)
        nqe.v = (double) phase[nqe.p] - TWOPI * int((wholepis+1)/2);
    else if (wholepis<=-1)
        nqe.v = (double) phase[nqe.p] + TWOPI * int((1-wholepis)/2);
    else
    	nqe.v = phase[nqe.p];
    
    unwrapped[nqe.p] = nqe.v;
    flag[nqe.p] = true;
    
    if (Push(primaryqueuenum,&nqe))
        std::cerr << "Out of memory!" << std::endl;

}

template <class T> void 
unwrap (const Matrix<size_t>& s, const size_t unb,
        const Matrix<std::complex<T> >& in, Matrix<T>& unwrapped) {
    
    Matrix<T> inamp = - abs(in);
    Matrix<T> wrapped = arg(in);
    unwrapped = Matrix<T> (in.Dims());

    size_t sze = numel(wrapped), i, sp, nb;
    container<size_t> dim = size(wrapped);

    // Minimum number of unwrapping bins is 2
    nb = (unb<2) ? 2 : unb;
    
    // Find min and max
    double minin, maxin, diff;
    minin    = min(inamp);
    maxin    = max(inamp);
    diff     = 1.00001 * (maxin - minin);
    
    sp = sub2ind(wrapped,s.Container());
    
    flag.resize(sze);
    
    for (i = 0; i < sze; i++)
        flag[i] = false;
    
    for (i = 0; i < nb; i++)
        InitQueue (i, sizeof(QUEUEENTRY));

    QUEUEENTRY qe;
    
    qe.p  = sp;
    qe.x  = s[0];
    qe.y  = s[1];
    qe.z  = s[2];
    unwrapped[85]  = wrapped[85];
    unwrapped[sp]  = qe.v = wrapped[sp];
    flag[sp] = true;
    
    // push s
    Push(0, &qe);
    
    // First, work out pole field threshold that we're going to use
    container<double> polefieldthresholds (nb);
    int ind  = 0;

    for (i = 0; i < nb; ++i)
        polefieldthresholds[ind++] = minin + diff * i / (nb-1);
    
    for (i=0; i < nb; i++) {

        while (!Pop(i,&qe)) {

            if (inamp[qe.p] > polefieldthresholds[i]) { // too close to a scary pole, so just defer by pushing to other stack
                ind = i;
                while (inamp[qe.p] > polefieldthresholds[++ind]);
                Push(ind,&qe);  // just defer by pushing to relevant stack
            } else { 
                Check (i, &qe,  m_bsz,  0,  0,  1, dim, wrapped, unwrapped);
                Check (i, &qe, -m_bsz,  0,  0, -1, dim, wrapped, unwrapped);
                Check (i, &qe,  m_bsy,  0,  1,  0, dim, wrapped, unwrapped);
                Check (i, &qe, -m_bsy,  0, -1,  0, dim, wrapped, unwrapped);
                Check (i, &qe,  m_bsx,  1,  0,  0, dim, wrapped, unwrapped);
                Check (i, &qe, -m_bsx, -1,  0,  0, dim, wrapped, unwrapped);
            }
        }

        TerminateQueue(i);	// done with this Queue so free up memory
        
    }
    
}


template<class T> Matrix<T>
unwrap3d (const Matrix<T>& M, const Matrix<size_t>& seed, const size_t nub = DEFAULT_NB) {

    assert (is3d(M));
    assert (seed[0] < size(M,0) && seed[1] < size(M,1) && seed[2] < size(M,2));

    Matrix<T> ret (size(M));
    m_bsx = M.Dszs()[0]; m_bsy = M.Dszs()[1]; m_bsz = M.Dszs()[2];
    unwrap (seed, nub, complex2(ones<T>(size(M)),M), ret);

    return ret;

}


template <class T> Matrix<T>
unwrap3d (const Matrix<T>& M, const size_t nub = DEFAULT_NB) {

    Matrix<size_t> seed (3,1);
    for (size_t i = 0; i < 3; i++)
        seed[i] = floor (.5 * (float)size(M,i)) -1;
    return unwrap3d (complex2(ones<T>(size(M)),M), seed, nub);

}


template <class T> Matrix<T>
unwrap3d (const Matrix<std::complex<T> >& M, const Matrix<size_t>& seed, const size_t nub = DEFAULT_NB) {
    
    assert (is3d(M));
    assert (seed[0] < size(M,0) && seed[1] < size(M,1) && seed[2] < size(M,2));

    Matrix<T> ret (size(M));
    m_bsx = M.Dszs()[0]; m_bsy = M.Dszs()[1]; m_bsz = M.Dszs()[2];
    unwrap (seed, nub, M, ret);
    
    return ret;
    
}


template <class T> Matrix<T>
unwrap3d (const Matrix<std::complex<T> >& M, const size_t nub = DEFAULT_NB) {

    Matrix<size_t> seed (3,1);
    for (size_t i = 0; i < 3; i++)
        seed[i] = floor (.5 * (float)size(M,i)) -1;
    return unwrap3d (M, seed, nub);
   
}
