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

#ifndef __CS_XSENSE_HPP__
#define __CS_XSENSE_HPP__

#include "CSENSE.hpp"
#include "NCSENSE.hpp"
#include "CX.hpp"
#include "NLCG.hpp"
#include "LBFGS.hpp"
#include "SplitBregman.hpp"
#include "DWT.hpp"
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
template <class T> class CS_XSENSE : public FT<T> {


public: 

    typedef typename TypeTraits<T>::RT RT;

    CS_XSENSE () : ft(0), dwt(0), nlopt(0) {/*TODO: Default constructor*/}
    virtual ~CS_XSENSE () {
        if (ft)
            delete ft;
        if (nlopt)
            delete nlopt;
        if (dwt)
            delete dwt;
        for (size_t i = 0; i < tvt.size(); ++i)
            delete tvt[i];
    }

    CS_XSENSE (const Params& p) : ft(0), dwt(0), nlopt(0) {
        
        printf ("Intialising CS_XSENSE ...\n");
        std::string key;

        _tvw.resize(2);
        _tvw[0] = try_to_fetch<float> (p, "tvw1", 0.);
        _tvw[1] = try_to_fetch<float> (p, "tvw2", 0.);
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
            printf ("**ERROR - CS_XSENSE: Invalid or unspecified FT operator\n");
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
            printf ("**ERROR - CS_XSENSE: Invalid or unspecified optimisation algorithm\n");
            throw UNDEFINED_OPTIMISATION_ALGORITHM;
            break;
        }
        
        _csiter = try_to_fetch<int> (p, "csiter", 0);
        _wf     = try_to_fetch<int> (p, "wl_family", -1);
        _wm     = try_to_fetch<int> (p, "wl_member", 0);

        if (_wf < -1 || _wf > 5)
            _wf = -1;

        if (_wf>-1)
            dwt = new DWT<T> (dims[0], (wlfamily)_wf, _wm);
        else 
            dwt = new DWT<T> (16, (wlfamily)_wf, _wm);

        tvt.PushBack(new TVOP<T>());
        tvt.PushBack(new TVOP<T>());
        
        printf ("... done.\n\n");

    }


   inline CS_XSENSE (const CS_XSENSE& tocopy) {
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
    

    inline virtual RT obj (const Matrix<T>& x, const Matrix<T>& dx, const RT& t, RT& rmse) {
        RT obj = Obj (t);
        rmse = sqrt(obj/_ndnz);
        if (_tvw[0])
            obj += TV (t,0);
        if (_tvw[1])
            obj += TV (t,1);
        if (_xfmw)
            obj += XFM (x, dx, t);
        return obj;
    }
    

    inline virtual Matrix<T> df (const Matrix<T>& x) {
        wx = *dwt->*x;
        Matrix<T> g = dObj (x);
        if (_xfmw)
            g += dXFM (x);
        if (_tvw[0])
            g += dTV  (x,0);
        if (_tvw[1])
            g += dTV  (x,1);
        return g;
    }
    
    inline virtual void Update (const Matrix<T>& dx) {
        wdx =  *dwt->*dx;        
        ffdbx = *ft * wx;
        ffdbg = *ft * wdx;
        if (_tvw[0]) {
            ttdbx = *tvt[0] * wx;
            ttdbg = *tvt[0] * wdx;
        }
    }

	virtual std::ostream& Print (std::ostream& os) const {
		FT<T>::Print(os);
        os << "    Weights: TV("<< _tvw[0] <<") TV("<< _tvw[1] <<") XF("<< _xfmw <<") L1("<<_l1<<") Pnorm: "
           <<_pnorm<< std::endl;
        os << *ft << std::endl;
        if (dwt)
            os << *dwt << std::endl;
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
        
        if (_verbose)
            vc.push_back(im_dc);

        _ndnz = (RT)nnz(data);

        im_dc  = *dwt * im_dc;

        printf ("  Running %i %s iterations ... \n", _csiter+1, nlopt_names[_nlopt_type]); fflush(stdout);

        for (size_t i = 0; i < (size_t)_csiter; i++) {
            nlopt->Minimise ((Operator<T>*)this, im_dc);
            if (_verbose)
                vc.push_back(*dwt ->* im_dc);
        }

        if (_verbose) {
            size_t cpsz = numel(im_dc);
            im_dc = zeros<T> (size(im_dc,0), size(im_dc,1), (_dim == 3) ?
                                  size(im_dc,2) : 1, vc.size());
            for (size_t i = 0; i < vc.size(); i++)
                std::copy (&vc[i][0], &vc[i][0]+cpsz, &im_dc[i*cpsz]);
        } else
            im_dc = *dwt ->* im_dc;
        
        return im_dc;

    }

    virtual FT<T>* getFT () {
        return ft;
    }
    
private:
    
    inline RT Obj (const RT& t) {
        om = ffdbx;
        if (t > 0.0)
            om += t * ffdbg;
        om -= data;
        return real(om.dotc(om));
    }
    
    inline RT TV (const RT& t, size_t i) {
        RT o = 0.0;
        om = ttdbx;
        if (t > 0.0)
            om += t * ttdbg;
        om *= conj(om);
        om += _l1;
        om ^= 0.5*_pnorm;
        for (size_t i = 0; i < om.Size(); i++)
            o += real(om[i]);
        return _tvw[i] * o;
    }
    
    inline RT XFM (const Matrix<T>& x, const Matrix<T>& g, const RT& t) {
        RT o = 0.;
        om = x;
        if (t > 0.0)
            om += t * g;
        om *= conj(om);
        om += _l1;
        om ^= 0.5*_pnorm;
        for (size_t i = 0; i < om.Size(); i++)
            o += om[i].real();
        return _xfmw * o;
    }
    
    
    /**
     * @brief Compute gradient of the data consistency
     */
    inline Matrix<T> dObj (const Matrix<T>& x) {
        return 2.0 * (*dwt * (*ft ->* ((*ft * wx) - data)));
    }
    
    
    /**
     * @brief Compute gradient of L1-transform operator
     *
     * @param  x   X
     * @param  cgp CG parameters
     * @return     The gradient
     */
    inline Matrix<T> dXFM (const Matrix<T>& x) {
        return _xfmw * (x * ((x * conj(x) + _l1) ^ (0.5*_pnorm-1.0)));
    }
    
    
    /**
     * @brief Compute gradient of the total variation operator
     *
     * @param  x   Image space original
     * @param  wx  Image space perturbance
     * @param  cgp Parameters
     */
    inline Matrix<T> dTV (const Matrix<T>& x, const size_t& i) {
        Matrix<T> dx, g;
        dx = *tvt[i] * wx;
        g  = dx * conj(dx);
        g += _l1;
        g ^= 0.5*_pnorm-1.0;
        g *= dx;
        g *= _pnorm;
        g  = *dwt * (*(tvt[i])->*g);
        return (_tvw[i] * g);
    }
    

    Params p;
    FT<T>* ft;
    DWT<T>* dwt;
    Vector<TVOP<T>* > tvt;
    NonLinear<T>* nlopt;
    Vector<size_t> _image_size;
    RT _xfmw, _l1, _pnorm, _tv1w, _tv2w;
    Vector<RT> _tvw;
    mutable RT _ndnz;
    int _verbose, _ft_type, _csiter, _wf, _wm, _nlopt_type, _dim;
    Matrix<T> ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx, om;
    mutable Matrix<T> data;

    
};

#endif
