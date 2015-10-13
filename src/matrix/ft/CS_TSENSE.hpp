/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
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

#ifndef __CS_TSENSE_HPP__
#define __CS_TSENSE_HPP__

#include "CSENSE.hpp"
#include "NCSENSE.hpp"
#include "CX.hpp"
#include "NLCG.hpp"
#include "LBFGS.hpp"
#include "SplitBregman.hpp"
#include "TVOP.hpp"

#include "Workspace.hpp"

#include <numeric>
#include <thread>

enum CS_EXCEPTION {UNDEFINED_FT_OPERATOR, UNDEFINED_OPTIMISATION_ALGORITHM};

const char* nlopt_names[] = {"NLCG", "L-BFGS", "Split Bregman"};

    using namespace codeare::optimisation;

/**
 * @brief Compressed sensing on Cartesian and non-Cartesian SENSE<br/>
 *        According Lustig et al.
 *
 */
template <class T> class CS_TSENSE : public FT<T> {


public: 

    typedef typename TypeTraits<T>::RT RT;

    CS_TSENSE () : ft(0), nlopt(0), _ft_type(0), _nlopt_type(0), _dim(2),
    		_csiter(0), _verbose(0){/*TODO: Default constructor*/}
    virtual ~CS_TSENSE () {
        if (ft)
            delete ft;
        if (nlopt)
            delete nlopt;
        for (size_t i = 0; i < tvt.size(); ++i)
            delete tvt[i];
    }

    CS_TSENSE (const Params& p) : ft(0), nlopt(0) {
        
        printf ("Intialising CS_TSENSE ...\n");
        std::string key;

        _tvw.resize(2);
        _tvv.resize(2);
        _tvw[0] = try_to_fetch<float> (p, "tvw1", 0.);
        _tvv[0] = try_to_fetch<Vector<size_t> > (p, "tv1", Vector<size_t>());
        _tvw[1] = try_to_fetch<float> (p, "tvw2", 0.);
        _tvv[1] = try_to_fetch<Vector<size_t> > (p, "tv2", Vector<size_t>());
        _xfmw = try_to_fetch<float> (p, "xfmw", 0.);
        _l1 = try_to_fetch<float> (p, "l1", 0.);
        _pnorm = try_to_fetch<float> (p, "pnorm", 0.);
        _image_size = try_to_fetch<Vector<size_t> > (p, "imsz", _image_size);
        _dim = _image_size.size();
        Vector<size_t> dims = try_to_fetch<Vector<size_t> > (p, "dims", _image_size);
        
        _verbose = try_to_fetch<int> (p, "verbose", 0.);
        _nlopt_type = try_to_fetch<int> (p, "nlopt", 0.);

        _ft_type = try_to_fetch<int> (p, "ft", 4);

        switch (_ft_type)
        {
		case 0: // Regular cartesian 
			ft = (FT<T>*) new DFT<T> (p);
			break;
        case 1: // Cartesian SENSE
			ft = (FT<T>*) new CSENSE<T> (p);
            break;
        case 2: // Regualar non-Cartesian 
			ft = (FT<T>*) new NFFT<T> (p);
            break;
        case 3: // Regualar non-Cartesian 
			ft = (FT<T>*) new NCSENSE<T> (p);
            break;
		default:
            printf ("**ERROR - CS_TSENSE: Invalid or unspecified FT operator\n");
			throw UNDEFINED_FT_OPERATOR;
			break;
        }

        switch (_nlopt_type)
        {
        case 0: // NLCG (demo only very slow)
            nlopt = (NonLinear<T>*) new NLCG<T> (p);
            break;
        case 1: // L-BFGS
            nlopt = (NonLinear<T>*) new LBFGS<T> (p);
            break;
        case 2: // Split-Bregman
            nlopt = (NonLinear<T>*) new SplitBregman<T> (p);
            break;
        default:
            printf ("**ERROR - CS_TSENSE: Invalid or unspecified optimisation algorithm\n");
            throw UNDEFINED_OPTIMISATION_ALGORITHM;
            break;
        }
        
        _csiter = try_to_fetch<int> (p, "csiter", 0);

        tvt.PushBack(new TVOP<T>(_tvv[0]));
        tvt.PushBack(new TVOP<T>(_tvv[1]));
        
        printf ("... done.\n\n");

    }


   inline CS_TSENSE (const CS_TSENSE& tocopy) {
        *this = tocopy;
    }


    /**
	 * @brief      Assign k-space trajectory
	 * 
	 * @param  k   K-space trajectory
	 */
    
	inline void KSpace (const Matrix<RT>& k) NOEXCEPT {
        ft->KSpace(k);
	}
	
    
	/**
	 * @brief      Assign k-space weigths (jacobian of k in t) 
	 * 
	 * @param  w   Weights
	 */
	inline void Weights (const Matrix<RT>& w) NOEXCEPT {
        ft->Weights(w);
	}

    /**
	 * @brief   Set k-space mask
	 * @param   mask  k-space mask
	 */
	inline void Mask (const Matrix<RT>& mask) NOEXCEPT {
        ft->Mask(mask);
	}
    

	virtual std::ostream& Print (std::ostream& os) const {
		FT<T>::Print(os);
        os << "    Weights: TV("<< _tvw[0] <<") TV("<< _tvw[1] <<") XF("<< _xfmw <<") L1("<<_l1<<") Pnorm: "
           <<_pnorm<< std::endl;
        os << *ft << std::endl;
        if (_tvw[0])
            os << *tvt[0] << std::endl;
        if (_tvw[1])
            os << *tvt[1] << std::endl;
        os << *nlopt ;
		return os;
	}


    inline virtual Matrix<T> operator* (const Matrix<T>& image_space) const {
        return Trafo(image_space);
    }

    inline virtual Matrix<T> Trafo (const Matrix<T>& m) const {
        return *ft * m;
    }

    inline virtual Matrix<T> operator->* (const Matrix<T>& k_space) const {
        return Adjoint(k_space);
    }

    inline virtual Matrix<T> Adjoint (const Matrix<T>& m) const {
        data = m;

        Matrix<T> im_dc;
        std::vector< Matrix<T> > vc;
                
        im_dc  = data;
        if (_ft_type != 2 && _ft_type != 3)
            im_dc /= wspace.Get<RT>("pdf");
        im_dc  = *ft ->* im_dc;
        
        _ndnz = (RT)nnz(data);

        RT ma = max(abs(im_dc));
        _tvw[0] *= ma;
        _tvw[1] *= ma;

        printf ("  Running %i %s iterations ... \n", _csiter, nlopt_names[_nlopt_type]); fflush(stdout);

        for (size_t i = 0; i < (size_t)_csiter; i++)
            nlopt->Minimise ((Operator<T>*)this, im_dc);

        return im_dc;

    }

    virtual FT<T>* getFT () {
        return ft;
    }
    
    inline virtual RT obj (const Matrix<T>& x, const Matrix<T>& dx, const RT& t, RT& rmse) const {
        RT obj = Obj (x,dx,t), objtv1 = 0, objtv2 = 0;
        rmse = sqrt(obj/_ndnz);
#pragma omp parallel num_threads(2)
        {
            size_t tn = omp_get_thread_num();
            switch (tn)
        	{
        	case 0:
        		if (_tvw[tn])
        			objtv1 = TV (x,dx,t,tn);
        		break;
        	case 1:
        		if (_tvw[tn])
        			objtv2 = TV (x,dx,t,tn);
        		break;
        	default: break;
        	}
        }
        return obj + objtv1 + objtv2;
    }
    
    inline virtual Matrix<T> df (const Matrix<T>& x) {
        Matrix<T> g = dObj (x);
        for (size_t i = 0; i < _tvw.size(); ++i)
            if (_tvw[i])
                g += dTV  (x,i);
        return g;
    }
    
    inline virtual void Update (const Matrix<T>& dx) {
    }

private:
    
    inline RT Obj (const Matrix<T>& x, const Matrix<T>& dx, const RT& t) const {
        Matrix<T> w = *ft * (x + t*dx) - data; 
        return real(w.dotc(w));
    }
    
    inline RT TV  (const Matrix<T>& x, const Matrix<T>& dx, const RT& t, size_t i) const {
        RT o = 0.;
        Matrix<T> w = *tvt[i] * (x + t*dx);
        w = (w*conj(w)+_l1)^.5;
        for (size_t i = 0; i < w.Size(); i++)
            o += real(w[i]);
        return _tvw[i] * o;
    }
    
    /**
     * @brief Compute gradient of the data consistency
     */
    inline Matrix<T> dObj (const Matrix<T>& x) const {
        return 2.0 * (*ft ->* ((*ft * x) - data));
    }
    
    
    /**
     * @brief Compute gradient of the total variation operator
     *
     * @param  x   Image space original
     * @param  wx  Image space perturbance
     * @param  cgp Parameters
     */
    inline Matrix<T> dTV (const Matrix<T>& x, const size_t& i) const {
        Matrix<T> w = *tvt[i] * x;
        return _tvw[i] * (*tvt[i]) ->* (w * ((w*conj(w)+_l1)^(-.5)));
    }
    

    Params p;
    FT<T>* ft;
    Vector<TVOP<T>* > tvt;
    Vector<Vector<size_t> > _tvv;
    NonLinear<T>* nlopt;
    Vector<size_t> _image_size;
    RT _xfmw, _l1, _pnorm;
    mutable Vector<RT> _tvw;
    mutable RT _ndnz;
    int _verbose, _ft_type, _csiter, _nlopt_type, _dim;
    Matrix<T> ffdbx, ffdbg, wx, wdx;
    Vector<Matrix<T> > ttdbx, ttdbg;
    mutable Matrix<T> data;

    
};

#endif
