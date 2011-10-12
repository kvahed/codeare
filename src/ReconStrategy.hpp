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
		ReconStrategy   () {}
		

		/**
		 * @brief       Default destructor
		 */ 
		virtual
		~ReconStrategy  () {}
		
		
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
		}
		

		/**
		 * @brief       Free RAM
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Finalise        () {
			
			while (!m_cplx.empty()) {
				cout << "Clearing RAM of " <<  m_cplx.begin()->first.c_str() << endl;
				FreeCplx(m_cplx.begin()->first.c_str());
			}
			
			while (!m_real.empty()) {
				cout << "Clearing RAM of " <<  m_real.begin()->first.c_str() << endl;
				FreeReal(m_real.begin()->first.c_str());
			}
			
			while (!m_pixel.empty()) {
				cout << "Clearing RAM of " <<  m_pixel.begin()->first.c_str() << endl;
				FreePixel(m_pixel.begin()->first.c_str());
			}
			
			return RRSModule::OK;

		}
		
		
		/**
		 * @brief       Add a complex matrix to complex matrix container
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<cplx>&
		AddCplx         (const string name, Ptr< Matrix<cplx> > m) {
			
			assert (m_cplx.find (name) == m_cplx.end());
			m_cplx.insert (pair< string, Ptr < Matrix<cplx> > >(name, m));
			return *m_cplx[name];
			
		}
		
		
		/**
		 * @brief       Get reference to complex data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<cplx>&
		GetCplx         (const string name) {
			
			return *m_cplx[name];
			
		}
		
		
		/**
		 * @brief       Get reference to complex data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<double>&
		GetReal          (const string name) {
			
			return *m_real[name];
			
		}
		
		
		/**
		 * @brief       Get reference to complex data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<short>&
		GetPixel          (const string name) {
			
			return *m_pixel[name];
			
		}
		

		/**
		 * @brief       Remove a complex matrix from complex container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreeCplx        (const string name) {
			
			map<string, Ptr< Matrix<cplx> > >::iterator it = m_cplx.find (name);
			
			if (it == m_cplx.end())
				return false;
			
			delete it->second;
			m_cplx.erase(it);
			
			return true;
			
		}
		
		
		/**
		 * @brief       Add a real matrix to real matrix container
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Reference to 
		 */
		Matrix<double>&
		AddReal         (const string name, Ptr< Matrix<double> > m) {
			
			assert(m_real.find (name) == m_real.end());
			m_real.insert (pair< string, Ptr < Matrix<double> > >(name, m));
			return *m_real[name];
			
		}
		

		/**
		 * @brief       Remove a real matrix from real container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreeReal        (const string name) {
			
			map<string, Ptr< Matrix<double> > >::iterator it = m_real.find (name);

			if (it == m_real.end())
				return false;

			delete it->second;
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
		Matrix<short>&
		AddPixel        (const string name, Ptr< Matrix<short> > m) {

			assert (m_pixel.find (name) == m_pixel.end());
			m_pixel.insert (pair<string, Ptr< Matrix<short> > > (name, m));
			return *m_pixel[name];
			
		}


		/**
		 * @brief       Remove a pixel matrix from pixel container
		 *
		 * @param  name Name
		 * @return      Success
		 */
		bool 
		FreePixel        (const string name) {
			
			map<string, Ptr< Matrix<short> > >::iterator it = m_pixel.find (name);

			if (it == m_pixel.end())
				return false;

			delete it->second;
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
		GetCplx       (const string name, cplx_data& c)   {
			
			if (m_cplx.find (name) == m_cplx.end())
				return;
			
			Ptr< Matrix<cplx> > tmp = m_cplx[name];
			
			for (int j = 0; j < INVALID_DIM; j++) {
				c.dims[j] = tmp->Dim(j);
				c.res[j]  = tmp->Res(j);
			}
			
			c.vals.length(2 * tmp->Size()); 

			memcpy (&c.vals[0], &tmp->At(0), c.vals.length() * sizeof(float));

			FreeCplx(name);
			
		}
		
		
		/**
		 * @brief        Set data for recon (Remote access)
		 *
		 * @param name   Name
		 * @param c      Cplx data
		 */
		void 
		SetCplx       (const string name, const cplx_data& c)   {
			
			Ptr< Matrix<cplx> > tmp;
			
			if (m_cplx.find (name) == m_cplx.end())
				m_cplx.insert (pair<string, Ptr< Matrix<cplx> > > (name, tmp = NEW (Matrix<cplx>())));
			else
			    tmp = m_cplx[name];
			
			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = c.dims[i];
				tmp->Res(i) = c.res[i];
			}
			
			tmp->Reset ();
			
			memcpy (&tmp->At(0), &c.vals[0], c.vals.length() * sizeof(float));

		}
		

		/**
		 * @brief        Get data from recon (Local access)
		 *
		 * @param  name  Name
		 * @param  m     Cplx data storage 
		 */
		void
		GetCplx       (const string name, Matrix<cplx>& m)   {

			if (m_cplx.find (name) == m_cplx.end())
				return;
			
			m = *m_cplx[name];
			
		}
		
		
		/**
		 * @brief        Set data for recon (Local access)
		 *
		 * @param  name  Name
		 * @param m      Complex data
		 */
		void 
		SetCplx       (const string name, Matrix<cplx>& m)   {
			
			if (m_cplx.find (name) == m_cplx.end())
				m_cplx.insert (pair<string, Ptr< Matrix<cplx> > > (name, NEW (Matrix<cplx>())));
			
			*m_cplx[name] = m;
			
		}
		

		/**
		 * @brief        Get data from recon (Remote access)
		 *
		 * @param  name  Name
		 * @param  r     Real data storage
		 */
		void 
		GetReal        (const string name, real_data& r)   {
			
			if (m_real.find (name) == m_real.end())
				return;
			
			Ptr< Matrix<double> > tmp = m_real[name];
			
			for (int j = 0; j < INVALID_DIM; j++) {
				r.dims[j] = tmp->Dim(j);
				r.res[j]  = tmp->Res(j);
			}
			
			r.vals.length(tmp->Size()); 
			
			memcpy (&r.vals[0], &tmp->At(0), tmp->Size() * sizeof(double));

			FreeReal(name);
			
		}
		

		/**
		 * @brief         Set data for recon (Local access)
		 * 
		 * @param  name  Name
		 * @param  r     Real data
		 */
		void 
		SetReal        (const string name, const real_data& r)   {
			
			Ptr< Matrix<double> > tmp;

			if (m_real.find (name) == m_real.end())
				m_real.insert (pair<string, Ptr< Matrix<double> > > (name, tmp = NEW( Matrix<double>())));
			else
				tmp = m_real[name];
			
			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = r.dims[i];
				tmp->Res(i) = r.res[i];
			}
			
			tmp->Reset ();
			
			memcpy (&tmp->At(0), &r.vals[0], tmp->Size() * sizeof(double));

		}
		
		
		/**
		 * @brief         Get data from recon
		 *
		 * @param  name  Name
		 * @param  m      Real data
		 */
		void
		GetReal         (const string name, Matrix<double>& m)   {
			
			if (m_real.find (name) == m_real.end())
				return;
			
			m = *m_real[name];
			
		}
		
	

		/**
		 * @brief         Set data for recon
		 *
		 * @param  name  Name
		 * @param  m      Real data storage
		 */
		void 
		SetReal         (const string name, Matrix<double>& m)   {
			
			if (m_real.find (name) == m_real.end())
				m_real.insert (pair<string, Ptr< Matrix<double> > > (name, NEW( Matrix<double>())));
			
			*m_real[name] = m;
			
		}
		

		/**
		 * @brief         Get data from recon
		 *
		 * @param  name  Name
		 * @param  p     Pixel data storage
		 */
		void 
		GetPixel          (const string name, pixel_data& p)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				return;
			
			Ptr< Matrix<short> > tmp = m_pixel[name];

			for (int j = 0; j < INVALID_DIM; j++) {
				p.dims[j] = tmp->Dim(j);
				p.res[j]  = tmp->Res(j);
			}
			
			p.vals.length(tmp->Size()); 
			
			memcpy (&p.vals[0], &tmp->At(0), tmp->Size() * sizeof(short));

			FreePixel(name);
			
		}
		
		
		/**
		 * @brief        Set data for recon
		 *
		 * @param  name  Name
		 * @param  p     Pixel data
		 */
		void 
		SetPixel         (const string name, const pixel_data& p)   {
			
			Ptr< Matrix<short> > tmp;

			if (m_pixel.find (name) == m_pixel.end())
				m_pixel.insert (pair<string, Ptr< Matrix<short> > > (name, tmp = NEW (Matrix<short>())));
			else
				tmp = m_pixel[name];

			for (int i = 0; i < INVALID_DIM; i++) {
				tmp->Dim(i) = p.dims[i];
				tmp->Res(i) = p.res[i];
			}
			
			tmp->Reset ();
			
			memcpy (&tmp->At(0), &p.vals[0], tmp->Size() * sizeof(short));

		}

		
		/**
		 * @brief        Get data from recon
		 *
		 * @param  name  Name
		 * @param  m     Pixel data storage
		 */
		void
		GetPixel         (const string name, Matrix<short>& m)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				return;
			
			 m = *m_pixel[name];
			
		}

		
		/**
		 * @brief        Set data for recon
		 *
		 * @param  name  Name
		 * @param  m     Pixel data
		 */
		void 
		SetPixel         (const string name, Matrix<short>& m)   {
			
			if (m_pixel.find (name) == m_pixel.end())
				m_pixel.insert (pair<string, Ptr< Matrix<short> > > (name, NEW (Matrix<short>())));
			
			*m_pixel[name] = m;
			
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
		
		map < string, Ptr< Matrix<cplx>   > > m_cplx;   /*!< @brief Complex data repository   */
		map < string, Ptr< Matrix<double> > > m_real;   /*!< @brief Real data repository      */
		map < string, Ptr< Matrix<short>  > > m_pixel;  /*!< @brief Integer data respository  */
		
		string                          m_name;         /*!< @brief Name                      */
		bool                            m_initialised;  /*!< @brief Reco is initialised       */
		
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


