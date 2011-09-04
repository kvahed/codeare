/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#ifndef __RECON_STRATEGY_HPP__
#define __RECON_STRATEGY_HPP__

#include "Matrix.hpp"
#include "Configurable.hpp"

#include "DllExport.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

#include <cstdlib>
#include <complex>
#include <stdint.h>

#include "cycle.h"            // FFTW cycle implementation


using namespace RRSModule;
using namespace std;

namespace RRServer {

	/**
	 * @brief Strategy for reconstruction strategies
	 *        Derive hereof to expand the reconstruction toolbox
	 *
	 */
	class ReconStrategy : public Configurable {
		
		
	public:
		

		/**
		 * @brief       Default constructor
		 */ 
		ReconStrategy   () {};
		

		/**
		 * @brief       Default destructor
		 */ 
		virtual
		~ReconStrategy  () {

			Finalise();

		};
		
		
		/**
		 * @brief       Mandatory implementation of actual data procession
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Process         () = 0;
		

		/**
		 * @brief       Mandatory implementation of initialiser
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Init            () = 0;
		

		/**
		 * @brief       Additional preparation?
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Prepare         () { 
			return RRSModule::OK; 
		};
		

		/**
		 * @brief       Free RAM
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Finalise        () {

			while (!m_cplx.empty()) {
				cout << "Clearing RAM of " <<  m_cplx.begin()->first.c_str() << endl;
				m_cplx.erase(m_cplx.begin());
			}

			while (!m_real.empty()) {
				cout << "Clearing RAM of " <<  m_real.begin()->first.c_str() << endl;
				m_real.erase(m_real.begin());
			}

			while (!m_pixel.empty()) {
				cout << "Clearing RAM of " <<  m_pixel.begin()->first.c_str() << endl;
				m_pixel.erase(m_pixel.begin());
			}


		};
		
		
		/**
		 * @brief       Add a complex matrix to complex matrix container
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		bool
		AddCplx         (const string name, Matrix<cplx>* m) {

			if (m_cplx.find (name) == m_cplx.end())
				m_cplx.insert (pair<string, Matrix<cplx>*> (name, m));
			else 
				return false;

			return true;
			
		}


		/**
		 * @brief       Remove a complex matrix from complex container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreeCplx        (const string name) {
			
			map<string,Matrix<cplx>*>::iterator it = m_cplx.find (name);

			if (it == m_cplx.end())
				return false;

			m_cplx.erase(it);

			return true;
			
		}
		

		/**
		 * @brief       Add a real matrix to real matrix container
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		bool
		AddReal         (const string name, Matrix<double>* m) {

			if (m_real.find (name) == m_real.end())
				m_real.insert (pair<string, Matrix<double>*> (name, m));
			else 
				return false;

			return true;
			
		}


		/**
		 * @brief       Remove a real matrix from real container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreeReal        (const string name) {
			
			map<string,Matrix<double>*>::iterator it = m_real.find (name);

			if (it == m_real.end())
				return false;

			m_real.erase(it);

			return true;
			
		}
		

		/**
		 * @brief       Add a pixel matrix to pixel matrix container
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		bool
		AddPixel        (const string name, Matrix<short>* m) {

			if (m_pixel.find (name) == m_pixel.end())
				m_pixel.insert (pair<string, Matrix<short>*> (name, m));
			else 
				return false;

			return true;
			
		}


		/**
		 * @brief       Remove a pixel matrix from pixel container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreePixel        (const string name) {
			
			map<string,Matrix<short>*>::iterator it = m_pixel.find (name);

			if (it == m_pixel.end())
				return false;

			m_pixel.erase(it);

			return true;
			
		}
		

		/**
		 * @brief        Get data from recon (Remote access)
		 *
		 * @param  name  Name
		 * @param  c     Raw data storage 
		 */
		void
		GetCplx       (const string name, cplx_data* c)   {
			
			if (m_cplx.find (name) == m_cplx.end())
				return;
			
			Matrix<cplx>* tmp = m_cplx[name];
			
			for (int j = 0; j < INVALID_DIM; j++) {
				c->dims[j] = tmp->Dim(j);
				c->res[j]  = tmp->Res(j);
			}
			
			c->dreal.length(tmp->Size()); 
			c->dimag.length(tmp->Size());
			
			for (int i = 0; i < tmp->Size(); i++) {
				c->dreal[i] = tmp->At(i).real();
				c->dimag[i] = tmp->At(i).imag(); 
			}
			
			tmp->Clear();
			
		}
		
		
		/**
		 * @brief        Set data for recon (Remote access)
		 *
		 * @param name   Name
		 * @param c      Cplx data
		 */
		void 
		SetCplx       (const string name, const cplx_data* c)   {

			Matrix<cplx>* tmp;
				
			if (m_cplx.find (name) == m_cplx.end())
				m_cplx.insert (pair<string, Matrix<cplx>*> (name, tmp = new Matrix<cplx>()));
			else
			    tmp = m_cplx[name];
			
			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = c->dims[i];
				tmp->Res(i) = c->res[i];
			}
			
			tmp->Reset ();
			
			for (int j = 0; j < tmp->Size(); j++)
				tmp->At(j) =  complex<float> (c->dreal[j], c->dimag[j]);
			
		}
		

		/**
		 * @brief        Get data from recon (Local access)
		 *
		 * @param  name  Name
		 * @param  m     Cplx data storage 
		 */
		void
		GetCplx       (const string name, Matrix<cplx>* m)   {

			if (m_cplx.find (name) == m_cplx.end())
				return;
			
			m = m_cplx[name];
			
		}
		
		
		/**
		 * @brief        Set data for recon (Local access)
		 *
		 * @param  name  Name
		 * @param m      Complex data
		 */
		void 
		SetCplx       (const string name, Matrix<cplx>* m)   {
			
			if (m_cplx.find (name) == m_cplx.end())
				m_cplx.insert (pair<string, Matrix<cplx>*> (name, new Matrix<cplx>()));
			
			m_cplx[name] = m;
			
		}
		

		/**
		 * @brief        Get data from recon (Remote access)
		 *
		 * @param  name  Name
		 * @param  r     Real data storage
		 */
		void 
		GetReal        (const string name, real_data* r)   {
			
			if (m_real.find (name) == m_real.end())
				return;
			
			Matrix<double>* tmp = m_real[name];
			
			for (int j = 0; j < INVALID_DIM; j++) {
				r->dims[j] = tmp->Dim(j);
				r->res[j]  = tmp->Res(j);
			}
			
			r->vals.length(tmp->Size()); 
			
			for (int i = 0; i < tmp->Size(); i++)
				r->vals[i] = tmp->At(i);
			
			tmp->Clear();
			
		}
		

		/**
		 * @brief         Set data for recon (Local access)
		 * 
		 * @param  name  Name
		 * @param  r     Real data
		 */
		void 
		SetReal        (const string name, const real_data* r)   {
			
			Matrix<double>* tmp;

			if (m_real.find (name) == m_real.end())
				m_real.insert (pair<string, Matrix<double>*> (name, tmp = new Matrix<double>()));
			else
				tmp = m_real[name];
			
			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = r->dims[i];
				tmp->Res(i) = r->res[i];
			}
			
			tmp->Reset ();
			
			for (int j = 0; j < tmp->Size(); j++)
				tmp->At(j) =  r->vals[j];
			
		}
		
		
		/**
		 * @brief         Get data from recon
		 *
		 * @param  name  Name
		 * @param  m      Real data
		 */
		void
		GetReal         (const string name, Matrix<double>* m)   {
			
			if (m_real.find (name) == m_real.end())
				return;
			
			m = m_real[name];
			
		}
		
	

		/**
		 * @brief         Set data for recon
		 *
		 * @param  name  Name
		 * @param  m      Real data storage
		 */
		void 
		SetReal         (const string name, Matrix<double>* m)   {
			
			if (m_real.find (name) == m_real.end())
				m_real.insert (pair<string, Matrix<double>*> (name, new Matrix<double>()));
			
			m_real[name] = m;
			
		}
		

		/**
		 * @brief         Get data from recon
		 *
		 * @param  name  Name
		 * @param  p     Pixel data storage
		 */
		void 
		GetPixel          (const string name, pixel_data* p)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				return;
			
			Matrix<short>* tmp = m_pixel[name];

			for (int j = 0; j < INVALID_DIM; j++) {
				p->dims[j] = tmp->Dim(j);
				p->res[j]  = tmp->Res(j);
			}
			
			p->vals.length(tmp->Size()); 
			
			for (int i = 0; i < tmp->Size(); i++)
				p->vals[i] = tmp->At(i);
			
			tmp->Clear();
			
		}
		
		
		/**
		 * @brief        Set data for recon
		 *
		 * @param  name  Name
		 * @param  p     Pixel data
		 */
		void 
		SetPixel         (const string name, const pixel_data* p)   {
			
			Matrix<short>* tmp;

			if (m_pixel.find (name) == m_pixel.end())
				m_pixel.insert (pair<string, Matrix<short>*> (name, tmp = new Matrix<short>()));
			else
				tmp = m_pixel[name];

			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = p->dims[i];
				tmp->Res(i) = p->res[i];
			}
			
			tmp->Reset ();
			
			for (int j = 0; j < tmp->Size(); j++) 
				tmp->At(j) =  p->vals[j];
			
		}

		
		/**
		 * @brief        Get data from recon
		 *
		 * @param  name  Name
		 * @param  m     Pixel data storage
		 */
		void
		GetPixel         (const string name, Matrix<short>* m)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				return;
			
			 m = m_pixel[name];
			
		}

		
		/**
		 * @brief        Set data for recon
		 *
		 * @param  name  Name
		 * @param  m     Pixel data
		 */
		void 
		SetPixel         (const string name, Matrix<short>* m)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				m_pixel.insert (pair<string, Matrix<short>*> (name, new Matrix<short>()));
			
			m_pixel[name] = m;
			
		}
	
		
		/**
		 * @brief          Attach a name to the algorithm
		 *
		 * @param  name    Name
		 */
		void Name (const char* name) { 

			m_name = string(name);
		
		}


		/**
		 * @brief          Get given name
		 *
		 * @return         Name
		 */
		const char* Name () {
			
			return m_name.c_str();
		
		}


	
	protected:
		
		map < string, Matrix< complex<float> >* >   m_cplx;
		map < string, Matrix<double>* >                  m_real;
		map < string, Matrix<short>* >                   m_pixel;
		
		string     m_name;        /*!< Name                               */
		
		bool            m_initialised; /*!< Reco is initialised                */
		
	};
	

}
#endif /* __RECON_STRATEGY_HPP__ */
	


/**
 * @brief              Dynamic constructor
 */
typedef RRServer::ReconStrategy* create_t  ();


/**
 *@brief               Dynamic destructor
 */
typedef void           destroy_t (RRServer::ReconStrategy*);


