/*
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

#include "IOFile.hpp"

#include "SiemensSoda.hpp"
#include "mdhVB15.hpp"

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
#include <stdlib.h>

#ifndef PARC_MODULE_NAME
enum IceDim {
    COL, LIN, CHA, SET, ECO, PHS, REP, SEG, PAR, SLC, IDA, IDB, IDC, IDD, IDE, AVE, MAX_ICE_DIM
};
#else
#define MAX_ICE_DIM 16
#endif


namespace codeare {
namespace matrix  {
namespace io      {

	enum EIM_BIT  {

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

	const size_t SYNC_HEADER_SIZE = 64;

	static bool bit_set (long eim, const EIM_BIT eb) {
		return (eim & (long)pow((float)2.,(float)eb));
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
	 * @brief  Data storage
	 */
	template<class T>
	struct Data {
		std::vector<size_t> dims;
		std::vector<size_t> idx;
		std::vector<float> ress;
		std::valarray<T>   data;

		/**
		 * @brief Constructor
		 */
		Data () {
			dims.resize(MAX_ICE_DIM,(size_t)1);
			ress.resize(MAX_ICE_DIM,(float)1);
			idx.resize(MAX_ICE_DIM,(size_t)1);
		}

		/**
		 * @brief Allocate RAM
		 */
		inline void
		Allocate (const bool verbose = false) {
			idx [0] = 1;
			for (size_t i = 1; i < MAX_ICE_DIM; i++)
				idx[i] = idx[i-1]*dims[i-1];
			printf ("  Data (dims: ");
			for (size_t i = 0; i < MAX_ICE_DIM; i++) {
				printf ("" JL_SIZE_T_SPECIFIER "", dims[i]);
				if (i < AVE)
					printf (" ");
				else
					printf (") (memory: %.1f MB)\n", 8.0*idx[AVE]/1024.0/1024.0);
			}
			data.resize(idx[AVE]);
		}

		/**
		 * @brief Clear RAM
		 */
		inline void
		Clear() {

			dims.resize(MAX_ICE_DIM,(size_t)1);
			ress.resize(MAX_ICE_DIM,(float)1);
			idx.resize(MAX_ICE_DIM,(size_t)1);
			data.resize(1, T(0));

		}

	};


	/**
	 * @brief    Reader for Syngo MR measurement file.
	 *           The reader provides a raw algorithm to transform the XProtocol to XML.
	 */
	class SyngoFile : public IOFile {

	public:

		typedef Data< std::complex<float> > ComplexData ;
		typedef Data< float > RealData ;

		SyngoFile (const std::string& fname, const IOMode& mode = READ,
				Params params = Params(), const bool& verbose = true) :
			IOFile(fname, mode, params, verbose), m_acqend(false),
			m_buffer(0), m_cur_sync(0), m_file(0), m_fod(0),
			m_meas_col_bytes(0), m_pos(0), m_size(0) {

			m_file   = ::fopen (fname.c_str(), "rb");

			if (m_file == NULL) {
				printf ("Error opening file '%s'\n", fname.c_str());
				m_status = FILE_OPEN_FAILED;
			}

			m_initialised = true;
			m_meas_dims = std::vector<size_t>(MAX_ICE_DIM,(size_t)1);

			dnames["NImageCols"] =  0; dnames["NLinMeas"] =  1; dnames["NSlcMeas"] =  2; dnames["NParMeas"] =  3;
			dnames[  "NEcoMeas"] =  4; dnames["NPhsMeas"] =  5; dnames["NRepMeas"] =  6; dnames["NSetMeas"] =  7;
			dnames[  "NSegMeas"] =  8; dnames[  "RawCha"] =  9; dnames["NIdaMeas"] = 10; dnames["NIdbMeas"] = 11;
			dnames[  "NIdcMeas"] = 12; dnames["NIddMeas"] = 13; dnames["NIdeMeas"] = 14; dnames["NAveMeas"] = 15;

		}


		error_code
		ReadFile (const bool verbose = true) {

			size_t read;

			m_buffer = (char*) malloc (sizeof(char)*m_size);

			if (m_buffer == NULL) {
				fputs ("Memory error", stderr);
				m_status = MEM_ALLOC_FAILED;
			}

			m_alloc = true;

			read = ::fread (m_buffer, sizeof(char), m_size, m_file);

			if (read != m_size) {
				fputs ("Reading error", stderr);
				m_status = FILE_READ_FAILED;
				return m_status;
			}

			if (verbose)
				printf ("  ... done reading %.1f MB\n\n", (float)(read)/1024.0/1024.0);

			m_pos = m_fod;

			return m_status;

		}


		void
		CalcSize (const bool verbose = true) {

			size_t read;
			unsigned fod;

			read = ::fread (&fod, sizeof(unsigned), 1, m_file);
			m_fod = (size_t) fod;
			if (read != 1) {
				fputs ("Reading error", stderr);
				m_status = FILE_READ_FAILED;
			}

			printf ("\n  Protocol:        %.1f KB\n", (-16.0+(float)m_fod)/1024.0);

			fseek (m_file,    -1, SEEK_END);
			m_size  = ftell(m_file);

			printf ("  Measurement:     %.1f MB\n", ((float)m_size-(float)m_fod)/1024.0/1024.0);
			fseek (m_file,     0, SEEK_SET);

		}


		error_code
		ParseProtocol (const bool verbose = true, const bool dry = false) {

			m_status = OK;

			return m_status;

		}


		error_code
		AllocateMatrices (const bool verbose = true) {

			m_status = OK;

			printf ("  Allocating data matrices ...\n");

			//std::copy(m_meas_dims.begin(), m_meas_dims.end(), std::ostream_iterator<int>(std::cout));

			float n = 8.0;
			for (size_t i = 0; i < m_meas_dims.size(); i++)
				n *= m_meas_dims[i];
			n /= pow(1024.0,3);

			printf ("\n%f\n", n);

			boost::any val = boost::make_shared<Matrix<cxfl> >(m_meas_dims);
			m_data.insert(std::pair<std::string, boost::any>("meas", val));
			printf ("  Data (dims: ");

			for (size_t i = 0; i < MAX_ICE_DIM; i++) {
				printf ("" JL_SIZE_T_SPECIFIER "", m_meas_dims[i]);
				if (i < AVE)
					printf (" ");
				else
					printf (")\n");// (memory: %.1f MB)\n", 8.0*idx[AVE]/1024.0/1024.0);
			}


			m_sync.Allocate();

			printf ("  ... done.\n\n");

			return m_status;

		}


		error_code
		ParseScan (const size_t& c, const bool dry = false) {

			m_status = OK;
			m_pos += m_meas_col_bytes;

			if (dry)
				return m_status;

            boost::shared_ptr<Matrix<cxfl> > m_meas =
					boost::any_cast<boost::shared_ptr<Matrix<cxfl> > >(m_data[std::string("meas")]);
			/*
			size_t dpos =
				m_cur_mdh.sLC.ushLine       * m_meas_dims[LIN] +
				m_cur_mdh.sLC.ushSlice      * m_meas_dims[SLC] +
				m_cur_mdh.sLC.ushPartition  * m_meas_dims[PAR] +
				m_cur_mdh.sLC.ushSeg        * m_meas_dims[SEG] +
				m_cur_mdh.sLC.ushEcho       * m_meas_dims[ECO] +
				m_cur_mdh.ushChannelId      * m_meas_dims[CHA] +
				m_cur_mdh.sLC.ushSet        * m_meas_dims[SET] +
				m_cur_mdh.sLC.ushPhase      * m_meas_dims[PHS] +
				m_cur_mdh.sLC.ushRepetition * m_meas_dims[REP] +
				m_cur_mdh.sLC.ushIda        * m_meas_dims[IDA] +
				m_cur_mdh.sLC.ushIdb        * m_meas_dims[IDB] +
				m_cur_mdh.sLC.ushIdc        * m_meas_dims[IDC] +
				m_cur_mdh.sLC.ushIdd        * m_meas_dims[IDD] +
				m_cur_mdh.sLC.ushIde        * m_meas_dims[IDE];

			memcpy (&m_meas->At(dpos), &m_buffer[m_pos], m_meas_col_bytes);
			*/

			printf ("(%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)\n",
					0,
					m_cur_mdh.sLC.ushLine,
					m_cur_mdh.sLC.ushSlice,
					m_cur_mdh.sLC.ushPartition,
					m_cur_mdh.sLC.ushSeg,
					m_cur_mdh.sLC.ushEcho,
					m_cur_mdh.ushChannelId,
					m_cur_mdh.sLC.ushSet,
					m_cur_mdh.sLC.ushPhase,
					m_cur_mdh.sLC.ushRepetition,
					m_cur_mdh.sLC.ushIda,
					m_cur_mdh.sLC.ushIdb,
					m_cur_mdh.sLC.ushIdc,
					m_cur_mdh.sLC.ushIdd,
					m_cur_mdh.sLC.ushIde);
			memcpy (&m_meas->At(0,
					m_cur_mdh.sLC.ushLine,
					m_cur_mdh.sLC.ushSlice,
					m_cur_mdh.sLC.ushPartition,
					m_cur_mdh.sLC.ushSeg,
					m_cur_mdh.sLC.ushEcho,
					m_cur_mdh.ushChannelId,
					m_cur_mdh.sLC.ushSet,
					m_cur_mdh.sLC.ushPhase,
					m_cur_mdh.sLC.ushRepetition,
					m_cur_mdh.sLC.ushIda,
					m_cur_mdh.sLC.ushIdb,
					m_cur_mdh.sLC.ushIdc,
					m_cur_mdh.sLC.ushIdd,
					m_cur_mdh.sLC.ushIde), &m_buffer[m_pos], m_meas_col_bytes);

			return m_status;

		}


		error_code
		ParseSyncData (const bool verbose = true, const bool dry = false) {

			int sync_size;

			m_status = OK;

			memcpy (&sync_size, &m_buffer[m_pos], sizeof(int));

			std::cout << sync_size << std::endl;

			m_pos += SYNC_HEADER_SIZE;

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


		error_code ParseMDH (const bool verbose = true) {

			m_status = OK;

			memcpy (&m_cur_mdh, &m_buffer[m_pos], sizeof(sMDH));

			if (bit_set(m_cur_mdh.aulEvalInfoMask[0], ACQEND)) {
				m_acqend = true;
				return m_status;
			}

			if (bit_set(m_cur_mdh.aulEvalInfoMask[0],SYNCDATA)) {
				m_sync.dims[1]++;
			} else if (m_meas_dims[0] == 1) {
				m_meas_dims[0] = (size_t)m_cur_mdh.ushSamplesInScan;
				m_meas_dims[6] = (size_t)m_cur_mdh.ushUsedChannels;
				m_meas_col_bytes = sizeof(std::complex<float>) * m_meas_dims[0];
			} else {
				if (m_meas_dims[LIN] < (size_t)m_cur_mdh.sLC.ushLine+1)
					m_meas_dims[LIN] = m_cur_mdh.sLC.ushLine+1;
				if (m_meas_dims[SLC] < (size_t)m_cur_mdh.sLC.ushSlice+1)
					m_meas_dims[SLC] = m_cur_mdh.sLC.ushSlice+1;
				if (m_meas_dims[PAR] < (size_t)m_cur_mdh.sLC.ushPartition+1)
					m_meas_dims[PAR] = m_cur_mdh.sLC.ushPartition+1;
				if (m_meas_dims[SEG] < (size_t)m_cur_mdh.sLC.ushSeg+1)
					m_meas_dims[SEG] = m_cur_mdh.sLC.ushSeg+1;
				if (m_meas_dims[ECO] < (size_t)m_cur_mdh.sLC.ushEcho+1)
					m_meas_dims[ECO] = m_cur_mdh.sLC.ushEcho+1;
				if (m_meas_dims[SET] < (size_t)m_cur_mdh.sLC.ushSet+1)
					m_meas_dims[SET] = m_cur_mdh.sLC.ushSet+1;
				if (m_meas_dims[PHS] <(size_t)m_cur_mdh.sLC.ushPhase+1)
					m_meas_dims[PHS] = m_cur_mdh.sLC.ushPhase+1;
				if (m_meas_dims[REP] <(size_t) m_cur_mdh.sLC.ushRepetition+1)
					m_meas_dims[REP] = m_cur_mdh.sLC.ushRepetition+1;
				if (m_meas_dims[IDA] < (size_t)m_cur_mdh.sLC.ushIda+1)
					m_meas_dims[IDA] = m_cur_mdh.sLC.ushIda+1;
				if (m_meas_dims[IDB] < (size_t)m_cur_mdh.sLC.ushIdb+1)
					m_meas_dims[IDB] = m_cur_mdh.sLC.ushIdb+1;
				if (m_meas_dims[IDC] < (size_t)m_cur_mdh.sLC.ushIdc+1)
					m_meas_dims[IDC] = m_cur_mdh.sLC.ushIdc+1;
				if (m_meas_dims[IDD] < (size_t)m_cur_mdh.sLC.ushIdd+1)
					m_meas_dims[IDD] = m_cur_mdh.sLC.ushIdd+1;
				if (m_meas_dims[IDE] < (size_t)m_cur_mdh.sLC.ushIde+1)
					m_meas_dims[IDE] = m_cur_mdh.sLC.ushIde+1;
			}

			m_pos += sizeof(sMDH);
			return m_status;

		}

		error_code ParseMeas (const bool verbose = true, const bool dry = false) {

			m_status   = OK;
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
					break;
					//ParseSyncData(verbose, dry);
				else
					ParseScan(verbose, dry);

				sc++;

			}

			printf ("  ... done.\n\n");

			return m_status;

		}


		virtual error_code
		CleanUp (const bool verbose = true) {

			m_status = OK;

			if (m_initialised) {

				if (m_alloc)
					free (m_buffer);

				::fclose(m_file);

				m_sync.Clear();
				m_noise.Clear();
				m_acs.Clear();

				m_initialised = false;

			}

			return m_status;

		}


		error_code
		ExtractProtocol (const bool verbose = true) {

			m_status = OK;
			m_xprot.append(&m_buffer[19],m_fod);
			return m_status;

		}


		error_code
		Read (const alloc_mode amode = SCAN_MEAS) {

			m_status = OK;

			/* Calculate data size */
			CalcSize (this->m_verb);

			/* Read data into buffer */
			ReadFile (this->m_verb);

			/* Exract protocol */
			ExtractProtocol();

			/* Access matrix sizes */
			if (amode == SCAN_MEAS)
				ParseMeas(this->m_verb, true);
			else
				ParseProtocol(this->m_verb);

			/* Allocate matrices */
			AllocateMatrices(this->m_verb);

			/* Parse data into matrices */
			ParseMeas(this->m_verb);

			/* Free buffer */
			CleanUp();

			return m_status;

		}


		template<class T> Matrix<T>
		Read (const std::string& name) const {
			Matrix<T> M;
			return M;
		}

		template<class T> bool
		Write (const Matrix<T>& M, const std::string& uri) {
			return false;
		}

		template<class T> Matrix<T>
		Read (const TiXmlElement* txe) const {
			Matrix<T> M;
			return M;
		}

		template<class T> bool
		Write (const Matrix<T>& M, const TiXmlElement* txe) {
			return false;
		}


		error_code
		Read (const std::string& dname = "") {
			return Read (SCAN_MEAS);
		}


		error_code
		Write (const std::string& dname = "") {
			return UNIMPLEMENTED_METHOD;
		}

		/**
		 * @brief Meant as a statuc mechanism
		 */
		~SyngoFile () {
			CleanUp();
		}

	private:
		

		bool        m_acqend;

		char*       m_buffer;
		
		ComplexData m_acs;
		ComplexData m_noise;

		RealData    m_sync;

		std::string m_xprot;

		FILE*       m_file;

		std::map < std::string, size_t > dnames;

		size_t      m_fod;
		size_t      m_size;
		size_t      m_pos;
		size_t      m_cur_sync;
		size_t      m_meas_col_bytes;

		sMDH        m_cur_mdh;

		std::vector<size_t> m_meas_dims;



		/**
		 * @brief Meant as a static mechanism
		 */
		SyngoFile  () :
			m_acqend(false),
			m_buffer(0),
			m_cur_sync(0),
			m_file(0),
			m_fod(0),
			m_size(0),
			m_pos(0),
			m_meas_col_bytes(0) {}


	};


} // namespace io
} // namespace matrix
} // namespace codeare


