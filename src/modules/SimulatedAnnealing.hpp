/*
 *  This file is part of codeare 
 *
 *  Copyright (C) 2007-2010 Kaveh Vahedipour
 *                          Forschungszentrum Juelich, Germany
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

#ifndef __SIMULATED_ANNEALING_HPP__
#define __SIMULATED_ANNEALING_HPP__

#include <gsl/gsl_rng.h>
#include <time.h>


/**
 * @brief  Travelling k-space-man with simulated annealing.
 */
class SimulatedAnnealing {

public:


	/**
	 * @brief           Construct
	 *
	 * @param  k        K-Space points
	 * @param  maxit    Max search iterations
	 * @param  st       Start temperature
	 * @param  ft       Final temperature
	 * @param  cr       Cooling rate 
	 * @param  verb     Verbosity
	 * @param  accworse Accept worse solutions with a certain probability
	 */
	SimulatedAnnealing (const Matrix<double>& k, const size_t& maxit = 1000, const double& st = 300.0, 
						const double& ft = 0.0, const double& cr = 0.95, const bool& verb = true, 
						const bool& accworse = false) {
		
		m_maxit = maxit;
		m_st    = st;
		m_ft    = ft;
		m_cr    = cr;
		m_verb  = verb;
		m_accworse = accworse;
		m_k     = k;

		gsl_rng_env_setup();
		m_rng = gsl_rng_alloc (gsl_rng_taus2);
		gsl_rng_set (m_rng, time(0));

		Initialise ();
		m_ct = m_st;

		printf ("Simulated annealing set up\n");

	}
	

	/**
	 * @brief           Clean up behind us
	 */
	~SimulatedAnnealing () {};
	

	/**
	 * @brief           Cool down
	 */
	void Cool () {

		printf ("  Cooling down ... \n");

		ticks cgstart = getticks();
		double len, rndn, dl, prob;
       	size_t i = 0, j = 0;

		while (i < m_maxit && m_ct > m_ft) {

			
			// Swap two sites
			SwapRange (m_sol, m_nsol); 
			len = Lengthiness (m_k, m_nsol);
			
			// Shorter than before: keep
			if (len <= m_len) {
				
				m_sol = m_nsol;
				m_len = len;
				
				if (m_verb) {
					if (j % 5 == 0 && j > 0)
						printf ("\n");
					printf ("    %04zu: %02.4f", i, m_len);
					j++;
				}

			} else if (m_accworse) {
				
				// accept worse solution with some probability
				dl   = len - m_len;
				rndn = gsl_rng_uniform (m_rng);
				prob = exp (- abs(dl) / m_ct);
				
				if (prob > rndn) {
					
					m_sol = m_nsol;
					m_len = len;
					
					if (m_verb) {
						if (j % 5 == 0 && j > 0)
							printf ("\n");
						printf ("    %04zu: %02.4f", i, m_len);
						j++;
					}
					
				}
				
		  	}
			
			if (m_bestlen > m_len) {
				m_best = m_sol;
				m_bestlen = m_len;
				m_bestn = i;
			}

			// decrease temperature
			m_ct -= m_cr;
			i++;
			
		}

		
		printf ("\n  ... done. Best trip length (%zu): %.4f, WTime: %.4f seconds.\n\n", 
				m_bestn, m_bestlen, elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());


	}



	/**
	 * @brief           Calculate the trajectory length 
	 * 
	 * @param  k        K-space points 
	 * @param  order    Order of positions
	 * @return          Trajectory length
	 */
	double Lengthiness (const Matrix<double>& k, const Matrix<size_t>& order) {
		
		double length = 0, x[3], y[3];
		size_t j, nr = numel(order), nx, ny;
		
		// calculate the total length of the traj
		
		for (j = 0; j < nr; j++) {
			
			nx = order[j];
			ny = order[(j == nr-1) ? 0 : j + 1];
			
			x[0] = k (nx,0);
			x[1] = k (nx,1);
			x[2] = k (nx,2);
			
			y[0] = k (ny,0);
			y[1] = k (ny,1);
			y[2] = k (ny,2);
			
			length += sqrt(pow(x[0]-y[0],2) + pow(x[1]-y[1],2) + pow(x[2]-y[2],2));
			
		}
		
		return length;		
		
    }	
	
	

	Matrix<double> GetSolution () const {

		Matrix<double> res = zeros<double> (size(m_k));
		size_t    nr       = numel (m_best);
		
		while (nr--)  {
			size_t npos = m_best[nr];
			res (nr,0) = m_k(npos,0);
			res (nr,1) = m_k(npos,1);
			res (nr,2) = m_k(npos,2);
		}
		
		return res;

	}

	
	/**
	 * @brief        Initialise with random 
	 */
	void  Initialise () {

		printf ("  Initialising ...\n");
		
		size_t nr = size(m_k,0);

		// create an array to hold initial solution and new solutions
		m_sol  = Matrix<size_t> (nr,1);
		m_nsol = Matrix<size_t> (nr,1);
		
		// and set it to a random permutation
		Permute (m_sol);
		
       	m_len = Lengthiness (m_k, m_sol);	
		m_bestlen = m_len;
		
		printf ("    Initial lenth: %.2f\n", m_len);
		
		if (m_verb) 
			MXDump (GetSolution(), "start.mat", "start");

		printf ("  ... done\n");
	} 


    /** 
	 * @brief       Random permutation of order 1..N
	 */
    void Permute (Matrix<size_t>& sol) {
		
		printf ("    Permuting ..."); fflush (stdout);

		size_t nr = numel (sol);
		
		// insert integers 1..N do not touch start / end (k-space center)
		for (size_t i = 1; i < nr-1; i++)
			sol[i]  = i+1;
		
		// shuffle
		for (size_t i = 1; i < nr-1; i++) {
			size_t r    = gsl_rng_uniform_int (m_rng, nr - 1) + 1;     // int between 1 and N
			size_t swap = sol[r];
			sol[r]      = sol[i];
			sol[i]      = swap;
		}
      
		printf (" done\n");
    }


	/**
	 * @brief        Swap a ranges order
	 */
    void SwapRange (const Matrix<size_t>& sol, Matrix<size_t>& nsol) {
		
		size_t i, x, y, pos, nr = numel(sol);
		
		// copy the current solution
		nsol = sol;
		
		// pick two places do not touch start 
		x = gsl_rng_uniform_int (m_rng, nr-1) + 1;
		y = gsl_rng_uniform_int (m_rng, nr-1) + 1;
	    
		// Swap
		if (y < x) {
			pos = x; x = y;	y = pos;
		} else if (x == y)
			return;
	
		// Create a temporary array that holds all the values between place1 and place2
		size_t len = y - x + 1;
		size_t cut [len];

		i = len;
		while (i--)
			nsol[x+len-1-i] = sol[x+i];
		
    }




private: 
	
	double         m_cr;      /**< @brief Cooling rate                 */
    double         m_st;      /**< @brief Start temperature            */
    double         m_ft;      /**< @brief Final temperature            */
    double         m_ct;      /**< @brief Current temperature          */

    double         m_len;     /**< @brief Current fitness              */
	double         m_bestlen;

    size_t         m_maxit;   /**< @brief Maximum number of iterations */
    size_t         m_nr;      /**< @brief Number of sites              */
	size_t         m_bestn;

    Matrix<size_t> m_sol;  /**< @brief Solution                     */
    Matrix<size_t> m_nsol; /**< @brief New solution                 */
    Matrix<size_t> m_best; /**< @brief Best solution                 */

	Matrix<double> m_k;  /**< @brief coordinates                  */

    bool           m_verb;  /**< @brief Verbose                      */
	bool           m_accworse; /**< @brief Randomly accept worse trip */

	gsl_rng* m_rng;
	
};


#endif
