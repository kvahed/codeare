/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifdef HAVE_CYCLE_H
   #include "cycle.h"
#endif

#include <valarray>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <complex>
#include <stdio.h>
#include <string.h>
#include <Matrix.hpp>


#include "mdhVB15.hpp"

namespace matrix {

	namespace io {

namespace SyngoMR {

enum EIM_BIT {
	
	ACQEND = 0,             //  0: End
	RTFEEDBACK,             //  1: Real-time feedback
	HPFEEDBACK,             //  2: HP feedback
	ONLINE,                 //  3: Handle online
	OFFLINE,                //  4: Handle offline
	SYNCDATA,               //  5: Synchroneous data
	
	LASTSCANINCONCAT = 8,   //  8: Last scan in concat

	RAWDATACORRECTION = 10, // 10: Correct raw data with correction factor 
	LASTSCANINMEAS,         // 11: Last scan in measurement
	SCANSCALEFACTOR,        // 12: Specific additional scale factor
	SECONDHADAMARPULSE,     // 13: 2nd RF excitation of HADAMAR
	REFPHASESTABSCAN,       // 14: Reference phase stabilisation scan
	PHASESTABSCAN,          // 15: Phase stabilisation scan
	D3FFT,                  // 16: Subject to 3D FFT
	SIGNREV,                // 17: Sign reversal
	PHASEFFT,               // 18: Perform PE FFT
	SWAPPED,                // 19: Swapped phase/readout direction
	POSTSHAREDLINE,         // 20: Shared line
	PHASCOR,                // 21: Phase correction data
	PATREFSCAN,             // 22: PAT reference data
	PATREFANDIMASCAN,       // 23: PAT reference and imaging data
	REFLECT,                // 24: Reflect line
	NOISEADJSCAN,           // 25: Noise adjust scan
	SHARENOW,               // 26: All lines are acquired from the actual and previous e.g. phases
	LASTMEASUREDLINE,       // 27: Last measured line of all succeeding e.g. phases
	FIRSTSCANINSLICE,       // 28: First scan in slice (needed for time stamps)
	LASTSCANINSLICE,        // 29: Last scan in slice  (      "       "       )
	TREFFECTIVEBEGIN,       // 30: Begin time stamp for TReff (triggered measurement)
	TREFFECTIVEEND          // 31: End time stamp for TReff (triggered measurement)
	
};


	static bool bit_set (long eim, const EIM_BIT eb) {

		return (eim & (long)pow(2.0,(float)eb));
		
	}


	enum meas_mode {
		TSE,
		ADJ_MDS_FRE,
		ADJ_FRE_MDS,
		ADJ_FRE,
		DRY,
		EPI,
		UTE
	};
	
	enum alloc_mode {
		SCAN_MEAS,
		SCAN_PROT
	};
	
	

	
	static const std::string root_path  = 
		"/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']";
	static const std::string rofov_path = 
		"/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamDouble[@name='RoFOV']/@value";
	static const std::string pefov_path = 
		"/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamDouble[@name='PeFOV']/@value";
	static const std::string slctn_path = 
		"/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamArray[@name='SliceThickness']/ParamDouble[1]/Precision[1]/@value";
	static const std::string rawch_path = 
		"/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online2']/ParamFunctor[@name='rawobjprovider']/ParamLong[@name='RawCha']/@value";

	/**
	 * @brief  Parameters
	 */
	struct Params {
		
		short int mode;
		
		bool      undo_os;
		bool      rep_wise;
		bool      average;
		bool      ext_ref_only;
		
		Params () {
			mode         = TSE;
			undo_os      = 1;
			rep_wise     = 0;
			average      = 1;
			ext_ref_only = 0;
		}
		
	};

	/**
	 * @brief  Grid size and resolution
	 */
	template<class T>
	struct Data {
		
		std::vector<size_t> dims;
		std::vector<size_t> idx;
		std::vector<float> ress;
		std::valarray<T>   data;

		Data () {
			dims.resize(INVALID_DIM,(size_t)1);
			ress.resize(INVALID_DIM,(float)1);
			idx.resize(INVALID_DIM,(size_t)1);
		}

		inline void 
		Allocate (const bool verbose = false) {
			idx [0] = 1;
			for (size_t i = 1; i < INVALID_DIM; i++)
				idx[i] = idx[i-1]*dims[i-1];
			printf ("  Data (dims: ");
			for (size_t i = 0; i < INVALID_DIM; i++) {
				printf ("%zu", dims[i]);
				if (i < AVE)
					printf (" ");
				else
					printf (") (memory: %.1f MB)\n", 8.0*idx[AVE]/1024.0/1024.0);
			}
			data.resize(idx[AVE]);
		}

		inline void 
		Clear() {
			
			dims.resize(INVALID_DIM,(size_t)1);
			ress.resize(INVALID_DIM,(float)1);
			idx.resize(INVALID_DIM,(size_t)1);
			data.resize(1, T(0));

		}
		
	};
		
	enum MFE {

		MF_OK,
		MF_FILE_OPEN_ERROR,
		MF_MEM_ALLOC_FAIL,
		MF_READ_FAIL

	};
	
		
	
	/**
	 * @brief    Reader for Syngo MR measurement file.
	 *           The reader provides a raw algorithm to transform the XProtocol to XML.
	 */
	class RawParser {

	public:
		
		typedef Data< std::complex<float> > ComplexData ;
		typedef Data< float > RealData ;

		RawParser (const std::string& fname, const bool verbose = true) :  
			m_acqend (false), m_initialised (false), m_allocated (false) {
			
			m_status  = MF_OK;
			
			m_fname   = fname;
			m_verbose = verbose;

			m_file    = fopen (fname.c_str(), "rb");
			
			if (m_file == NULL) {
				printf ("Error opening file '%s'\n", fname.c_str());
				m_status = MF_FILE_OPEN_ERROR;
			}

			m_initialised = true;
			
			dnames["NImageCols"] =  0; dnames["NLinMeas"] =  1; dnames["NSlcMeas"] =  2; dnames["NParMeas"] =  3;
			dnames[  "NEcoMeas"] =  4; dnames["NPhsMeas"] =  5; dnames["NRepMeas"] =  6; dnames["NSetMeas"] =  7;
			dnames[  "NSegMeas"] =  8; dnames[  "RawCha"] =  9; dnames["NIdaMeas"] = 10; dnames["NIdbMeas"] = 11;
			dnames[  "NIdcMeas"] = 12; dnames["NIddMeas"] = 13; dnames["NIdeMeas"] = 14; dnames["NAveMeas"] = 15; 

		}

		
		MFE ReadFile (const bool verbose = true) {

			size_t read;

			m_buffer = (char*) malloc (sizeof(char)*m_size);

			if (m_buffer == NULL) {
				fputs ("Memory error", stderr); 
				m_status = MF_MEM_ALLOC_FAIL;
			}
			
			m_allocated = true;

			read = fread (m_buffer, sizeof(char), m_size, m_file);

			if (read != m_size) {
				fputs ("Reading error", stderr); 
				m_status = MF_READ_FAIL;
				return m_status;
			}

			if (verbose)
				printf ("  ... done reading %.1f MB\n\n", (float)(read)/1024.0/1024.0);

			m_pos = m_fod;

			return m_status;
			
		}

		
		void CalcSize (const bool verbose = true) {

			size_t read;
			unsigned fod;

			read = fread (&fod, sizeof(unsigned), 1, m_file);
			m_fod = (size_t) fod;
			if (read != 1) {
				fputs ("Reading error", stderr); 
				m_status = MF_READ_FAIL;
			}

			printf ("\n  Protocol:        %.1f KB\n", (-16.0+(float)m_fod)/1024.0);

			fseek (m_file,    -1, SEEK_END);
			m_size  = ftell(m_file);

			printf ("  Measurement:     %.1f MB\n", ((float)m_size-(float)m_fod)/1024.0/1024.0);
			fseek (m_file,     0, SEEK_SET);

		}
		

		MFE ParseProtocol (const bool verbose = true, const bool dry = false) {
			
			m_status = MF_OK;
			
			return m_status;
			
		}
		

		MFE AllocateMatrices (const bool verbose = true) {
			
			m_status = MF_OK;

			printf ("  Allocating data matrices ...\n");

			m_meas.Allocate();
			m_sync.Allocate();

			printf ("  ... done.\n\n");

			return m_status;
			
		}
		

		MFE ParseScan (const size_t& c, const bool dry = false) {

			m_status = MF_OK;
			m_pos += m_meas_col_bytes;

			if (dry) 
				return m_status;

			size_t dpos =  
				m_cur_mdh.sLC.ushLine       * m_meas.dims[LIN] +
				m_cur_mdh.sLC.ushSlice      * m_meas.dims[SLC] +
				m_cur_mdh.sLC.ushPartition  * m_meas.dims[PAR] +
				m_cur_mdh.sLC.ushSeg        * m_meas.dims[SEG] +
				m_cur_mdh.sLC.ushEcho       * m_meas.dims[ECO] +
				m_cur_mdh.ushChannelId      * m_meas.dims[CHA] +
				m_cur_mdh.sLC.ushSet        * m_meas.dims[SET] +
				m_cur_mdh.sLC.ushPhase      * m_meas.dims[PHS] +
				m_cur_mdh.sLC.ushRepetition * m_meas.dims[REP] +
				m_cur_mdh.sLC.ushIda        * m_meas.dims[IDA] +
				m_cur_mdh.sLC.ushIdb        * m_meas.dims[IDB] +
				m_cur_mdh.sLC.ushIdc        * m_meas.dims[IDC] +
				m_cur_mdh.sLC.ushIdd        * m_meas.dims[IDD] +
				m_cur_mdh.sLC.ushIde        * m_meas.dims[IDE];
			
			memcpy (&m_meas.data[dpos], &m_buffer[m_pos], m_meas_col_bytes);
			
			return m_status;
					
		}
		

		MFE ParseSyncData (const bool verbose = true, const bool dry = false) {

			int sync_size;
			
			m_status = MF_OK;
			
			memcpy (&sync_size, &m_buffer[m_pos], sizeof(int));
			
			m_pos += 64;

			if (dry) {
				if (m_sync.dims[1] < (size_t)sync_size/4)
					m_sync.dims[0] = sync_size/4;
			} else {
				memcpy (&m_sync.data[0], &m_buffer[m_pos], sync_size);
				m_cur_sync ++;
			}

			m_pos += sync_size;

			return m_status;

		}


		MFE ParseMDH (const bool verbose = true) {

			m_status = MF_OK;
			
			memcpy (&m_cur_mdh, &m_buffer[m_pos], sizeof(sMDH));

			if (bit_set(m_cur_mdh.aulEvalInfoMask[0], ACQEND)) {
				m_acqend = true;
				return m_status;
			}

			if (bit_set(m_cur_mdh.aulEvalInfoMask[0],SYNCDATA)) {
				m_sync.dims[1]++;
			} else if (m_meas.dims[0] == 1) {
				m_meas.dims[0] = (size_t)m_cur_mdh.ushSamplesInScan;
				m_meas.dims[6] = (size_t)m_cur_mdh.ushUsedChannels;
				m_meas_col_bytes = sizeof(std::complex<float>) * m_meas.dims[0];
			} else {
				if (m_meas.dims[LIN] < (size_t)m_cur_mdh.sLC.ushLine+1)
					m_meas.dims[LIN] = m_cur_mdh.sLC.ushLine+1;
				if (m_meas.dims[SLC] < (size_t)m_cur_mdh.sLC.ushSlice+1)
					m_meas.dims[SLC] = m_cur_mdh.sLC.ushSlice+1;
				if (m_meas.dims[PAR] < (size_t)m_cur_mdh.sLC.ushPartition+1)
					m_meas.dims[PAR] = m_cur_mdh.sLC.ushPartition+1;
				if (m_meas.dims[SEG] < (size_t)m_cur_mdh.sLC.ushSeg+1)
					m_meas.dims[SEG] = m_cur_mdh.sLC.ushSeg+1;
				if (m_meas.dims[ECO] < (size_t)m_cur_mdh.sLC.ushEcho+1)
					m_meas.dims[ECO] = m_cur_mdh.sLC.ushEcho+1;
				if (m_meas.dims[SET] < (size_t)m_cur_mdh.sLC.ushSet+1)
					m_meas.dims[SET] = m_cur_mdh.sLC.ushSet+1;
				if (m_meas.dims[PHS] <(size_t)m_cur_mdh.sLC.ushPhase+1)
					m_meas.dims[PHS] = m_cur_mdh.sLC.ushPhase+1;
				if (m_meas.dims[REP] <(size_t) m_cur_mdh.sLC.ushRepetition+1)
					m_meas.dims[REP] = m_cur_mdh.sLC.ushRepetition+1;
				if (m_meas.dims[IDA] < (size_t)m_cur_mdh.sLC.ushIda+1)
					m_meas.dims[IDA] = m_cur_mdh.sLC.ushIda+1;
				if (m_meas.dims[IDB] < (size_t)m_cur_mdh.sLC.ushIdb+1)
					m_meas.dims[IDB] = m_cur_mdh.sLC.ushIdb+1;
				if (m_meas.dims[IDC] < (size_t)m_cur_mdh.sLC.ushIdc+1)
					m_meas.dims[IDC] = m_cur_mdh.sLC.ushIdc+1;
				if (m_meas.dims[IDD] < (size_t)m_cur_mdh.sLC.ushIdd+1)
					m_meas.dims[IDD] = m_cur_mdh.sLC.ushIdd+1;
				if (m_meas.dims[IDE] < (size_t)m_cur_mdh.sLC.ushIde+1)
					m_meas.dims[IDE] = m_cur_mdh.sLC.ushIde+1;
			}

			m_pos += sizeof(sMDH);
			return m_status;
			
		}

		MFE ParseMeas (const bool verbose = true, const bool dry = false) {

			m_status   = MF_OK;
			m_acqend   = false;
			m_pos      = m_fod;
			m_cur_sync = 0;

			printf ("  Parsing measurement ");
			if (dry)
				printf ("(dry run!)");
			printf (" ...\n");

			size_t sc = 0;

			while (true) {

				ParseMDH(verbose);
				
				if (m_acqend)
					break;

				if (bit_set(m_cur_mdh.aulEvalInfoMask[0], SYNCDATA))
					ParseSyncData(verbose, dry);
				else 
					ParseScan(verbose, dry);

				sc++;

			}

			printf ("  ... done.\n\n");

			return m_status;

		}


		MFE CleanUp (const bool verbose = true) {

			m_status = MF_OK;

			if (m_initialised) {
				
				if (m_allocated)
					free (m_buffer);

				fclose(m_file);

				m_meas.Clear();
				m_sync.Clear();
				m_noise.Clear();
				m_acs.Clear();

				m_initialised = false;

			}

			return m_status;

		}
		
		
		MFE ExtractProtocol (const bool verbose = true) {

			m_status = MF_OK;
			m_xprot.append(&m_buffer[19],m_fod);
			return m_status;
			
		}


		MFE Parse (const alloc_mode amode = SCAN_MEAS, const bool verbose = true) {
			
			m_status = MF_OK;
			
			/* Calculate data size */
			CalcSize (verbose);

			/* Read data into buffer */
			ReadFile (verbose);
			
			/* Exract protocol */
			ExtractProtocol();

			/* Access matrix sizes */
			if (amode == SCAN_MEAS)
				ParseMeas(verbose, true);
			else
				ParseProtocol(verbose);

			/* Allocate matrices */
			AllocateMatrices(verbose);

			/* Parse data into matrices */
			ParseMeas(verbose);

			/* Free buffer */
			CleanUp();

			return m_status;

		}


		MFE Status () const {

			return m_status;

		}

		
		/**
		 * @brief Meant as a static mechanism
		 */
		RawParser  () : 
			m_acqend (false),
			m_initialised (false),
			m_allocated (false) {}
		

		/**
		 * @brief Meant as a statuc mechanism
		 */
		~RawParser () {

			CleanUp();

		}


		Data<cxfl> 
		GetMeas () const {

			return m_meas;

		}


	private:
		

		bool        m_acqend;
		bool        m_initialised;
		bool        m_allocated;

		char*       m_buffer;
		
		ComplexData m_meas;
		RealData    m_sync;
		ComplexData m_acs;
		ComplexData m_noise;
		
		MFE         m_status;
		
		std::string m_xprot;

		FILE*       m_file;
		std::string m_fname;
		bool        m_verbose;
		
		std::map < std::string, size_t > dnames;
	
		size_t      m_fod;
		size_t      m_size;
		size_t      m_pos;
		size_t      m_cur_sync;
		size_t      m_meas_col_bytes;

		sMDH        m_cur_mdh;
	

	};


}

	}
}
