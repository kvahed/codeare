/*
 *  This file is part of rawpp.
 *
 *  rawpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  rawpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rawpp. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 17 Jul, 2014
 *      Author: Kaveh Vahedipour
 */

#ifndef _VD_FILE_HPP_
#define _VD_FILE_HPP_

#include "SyngoFile.hpp"
namespace codeare {
namespace matrix  {
namespace io      {
namespace VD {
class VDFile :
    public SyngoFile {
    
public:

#ifdef USE_MATLAB
    // Construct 
    VDFile (const std::string fname, int nlhs = 0, mxArray *lhs[] = 0) :
        SyngoFile(fname, nlhs, lhs), _id(0), _ndset(1), _meas_r(0), _meas_i(0), _sync_r(0),
		_nmeas(0), _nlines(0), _nsync(0), _meas(0), _digested(false) {
        _file.seekg(0);
        _file.read ((char*)&_id, sizeof(uint32_t));      // ID
        _file.read ((char*)&_ndset, sizeof(uint32_t));   // # data sets
        prtmsg ("     VD raid file (id:%d) contains %d data set(s):\n", _id, _ndset);
        _veh.resize(_ndset);
        for (size_t i = 0; i < _ndset; ++i) {
            _file.read ((char*)&_veh[i], ENTRY_HEADER_LEN);
            prtmsg ("       MID(%d) FID(%d) %s \n", _veh[i].MeasID, _veh[i].FieldID, _veh[i].ProtocolName);
        }
        _measdims.resize(16);       // ICE dimensions
        _measdims_sizes.resize(14); // Sizes except for COL/LIN
        _syncdims.resize(2,0);      // Syncdata dimension
    }
#else
    VDFile (const std::string& fname, const IOMode mode, const Params& params, const bool verbosity) :
    	SyngoFile(fname, mode, params, verbosity), _id(0), _ndset(1), _meas_r(0), _meas_i(0),
		_sync_r(0), _nmeas(0), _nlines(0), _nsync(0), _digested(false), _tend(0), _tstart(0),
		_cent_par(0), _cent_col(0), _cent_lin(0), _ta(0), _tr(0) {
    	std::cout << _protocol.Get<long>("XProtocol.ASCCONV.lTotalScanTimeSec") << std::endl;
        _file.seekg(0);
        _file.read ((char*)&_id, sizeof(uint32_t));      // ID
        _file.read ((char*)&_ndset, sizeof(uint32_t));   // # data sets
        prtmsg ("     VD raid file (id:%d) contains %d data set(s):\n", _id, _ndset);
        _veh.resize(_ndset);
        for (size_t i = 0; i < _ndset; ++i) {
            _file.read ((char*)&_veh[i], ENTRY_HEADER_LEN);
            prtmsg ("       MID(%d) FID(%d) %s \n", _veh[i].MeasID, _veh[i].FieldID, _veh[i].ProtocolName);
        }
        _measdims.resize(16);       // ICE dimensions
        _measdims_sizes.resize(14); // Sizes except for COL/LIN
        _syncdims.resize(2,0);      // Syncdata dimension
    }
#endif


    // Clean up
    virtual ~VDFile () {
        _meas_r = 0;
        _meas_i = 0;
        _sync_r = 0;
        _measdims.clear();
    }


    /**
     * @brief Digest ingredients
     *
     * @param  dry    Dry run
     */
    virtual void Digest () {
        
        uint32_t cur_pos;
        MeasHeader sh;
        SimpleTimer st;
        
        if (_meas_r == 0)
            prtmsg ("   Parsing ... \n");
        else 
            prtmsg ("   Reading ... \n");
        
        _file.seekg(_veh.back().MeasOffset);                // Go to last measurement
        _file.read((char*)&cur_pos, sizeof(uint32_t));
        _file.seekg (_veh.back().MeasOffset + cur_pos);     // Skip protocol
        _nmeas = 0;
        _nlines = 0;
        
        while (true) {
            _file.read((char*)&sh, MEAS_HEADER_LEN);
            if (bit_set(sh.aulEvalInfoMask[0], ACQEND) || sh.ushSamplesInScan == 0) { // ACQEND
                break;
            } else if (bit_set (sh.aulEvalInfoMask[0], SYNCDATA)) {
                size_t tmp1 = _file.tellg(), tmp2;
                ParseSync (sh);
                tmp2 = _file.tellg();
                tmp1 = ((tmp2-tmp1)%32);
                if (tmp1)
                    _file.seekg (tmp2 + 32-tmp1);
            } else if (bit_set(sh.aulEvalInfoMask[1], ONLINE)) { // CT_NORMALIZE
                    ParseCTNormalize (sh);
            } else if (bit_set (sh.aulEvalInfoMask[0], ONLINE)) { // Actual data
                if (_nmeas == 0) {
                    _measdims[0] = sh.ushSamplesInScan;
                    _measdims[1] = sh.ushUsedChannels;
                }
                ParseMeas (sh);
            } else {
                prtwrn ("Failed to understand data set. Skipping.");
                _measdims.resize(16,1);
                break;
            }
            _nlines++;
        }
        prtmsg ("     done - wtime %s", st.Format().c_str());

        //prtmsg ("     done - wtime (%d ms).\n", 0/*GetTimeMs64()-start*/);
        
        if (_meas_r == 0) {
            _measdims = _raise_one (_measdims); // Dimensions one higher that highest counter
            _ta = 2.5e-3*(_tend-_tstart);
            wspace.PSet("TA", _ta);
            _tr = _ta/_nmeas;                            
            wspace.PSet("TR", _tr);
            PrintParse(); 
            Allocate(); 
            Digest(); 
        }
        _digested = true;
    }
    
    void Read () {
        if (!_digested)
            this->Digest();
        wspace.Add<cxfl>("meas", _meas);
        wspace.Add<cxfl>("rtfb", _rtfb);
        wspace.Add<float>("sync", _sync);
    }

private:

    /**
     * @brief Parse one multichannel measurement line
     *
     * @param  mh    Measurement header
     */
    inline void ParseMeas (const MeasHeader& mh) {
        size_t pos (_file.tellg());

        ChannelHeader ch;
        if (!_meas_r) {
            _measdims = _max (_measdims, mh.sLC);  // Keep track of highest counter
            _file.seekg (pos + (mh.ushSamplesInScan*sizeof(std::complex<float>)
                             + CHANNEL_HEADER_LEN)*mh.ushUsedChannels);
            if (_nmeas == 0) {
                _tstart   = mh.ulTimeStamp;
                _cent_par = mh.ushKSpaceCentrePartitionNo;
                _cent_lin = mh.ushKSpaceCentreLineNo;
                _cent_col = mh.ushKSpaceCentreColumn;
            } else if (bit_set(mh.aulEvalInfoMask[0], LASTSCANINMEAS)) {
                _tend = mh.ulTimeStamp;
            }

        } else {
#ifndef USE_IN_MATLAB
            for (size_t i = 0; i < mh.ushUsedChannels; ++i) {
                _file.read((char*)&ch, CHANNEL_HEADER_LEN);
                _file.read ((char*)&_meas(0, i, mh.sLC[0], mh.sLC[1], mh.sLC[2],
                    mh.sLC[3], mh.sLC[4], mh.sLC[5], mh.sLC[6], mh.sLC[7],
					mh.sLC[8], mh.sLC[9], mh.sLC[10], mh.sLC[11], mh.sLC[12],
					mh.sLC[13]), mh.ushSamplesInScan*sizeof(std::complex<float>));
            }
#else
            std::vector<std::complex<float> > buf (mh.ushSamplesInScan*mh.ushUsedChannels);
            for (size_t i = 0; i < mh.ushUsedChannels; ++i) {
                _file.read((char*)&ch, CHANNEL_HEADER_LEN);
                _file.read((char*)&buf[i*mh.ushSamplesInScan],
                           mh.ushSamplesInScan*sizeof(std::complex<float>));
            }
            size_t offset = 0;
            for (size_t i = 0; i < 14; ++i)
                offset += mh.sLC[i]*_measdims_sizes[i];
            for (size_t j = 0; j < mh.ushUsedChannels; ++j) {
                for (size_t i = 0; i < mh.ushSamplesInScan; ++i) {
                    _meas_r[offset + j*mh.ushSamplesInScan + i] =
                        std::real(buf[j*mh.ushSamplesInScan + i]);
                    _meas_i[offset + j*mh.ushSamplesInScan + i] =
                        std::imag(buf[j*mh.ushSamplesInScan + i]);
                }
            }
#endif
        }
        _nmeas++;
    }


    inline void ParseCTNormalize (const MeasHeader& mh) {
        ChannelHeader ch;
        std::vector<std::complex<float> > buf (mh.ushSamplesInScan*mh.ushUsedChannels);
        for (size_t i = 0; i < mh.ushUsedChannels; ++i) {
            _file.read((char*)&ch, CHANNEL_HEADER_LEN);
            _file.read((char*)&buf[i*mh.ushSamplesInScan],
                       mh.ushSamplesInScan*sizeof(std::complex<float>));
        }
    }
        
    /**
     * @brief Parse SyncData
     */
    inline void ParseSync (const MeasHeader& mh) {
        size_t start (_file.tellg());
        if (!_sync_r) {
            if (_syncdims[0] == 0) {
                _file.read((char*)&_syncdims[0], sizeof(uint32_t));
                _syncdims[0] /= sizeof(float);
            }
            _file.seekg (start + SYNC_HEADER_SIZE - 4 + _syncdims[0]*sizeof(float));
        } else {
            _file.seekg (start + SYNC_HEADER_SIZE - 4);
            _file.read((char*)&_sync_r[_syncdims[0]*_syncdims[1]], _syncdims[0]*sizeof(float));
        }
        _syncdims[1]++;
    }
    
    
    /**
     * @brief Allocate RAM for data
     */
    void Allocate () {

        msize_t measdim[16] = {
            _measdims[ 0], _measdims[ 1], _measdims[ 2], _measdims[ 3],
			_measdims[ 4], _measdims[ 5], _measdims[ 6], _measdims[ 7],
			_measdims[ 8], _measdims[ 9], _measdims[10], _measdims[11],
			_measdims[12], _measdims[13], _measdims[14], _measdims[15]};
        msize_t syncdim[2] = {_syncdims[0], _syncdims[1]};
        
#ifdef USE_IN_MATLAB
        _lhs[0] = mxCreateNumericArray (16, measdim, mxSINGLE_CLASS, mxCOMPLEX);
        _measdims_sizes[0] = _measdims[ 0]*_measdims[ 1];
        for (size_t i = 1; i < 14; ++i)
            _measdims_sizes[i] = _measdims_sizes[i-1]*_measdims[i+1];
        _meas_r = (float*) mxGetData(_lhs[0]);
        _meas_i = (float*) mxGetImagData(_lhs[0]);
        _nmeas = 0;
        if (_nlhs > 1) {
            _lhs[1] = mxCreateNumericArray (2, syncdim, mxSINGLE_CLASS, mxREAL);
            _sync_r = (float*) mxGetData(_lhs[1]);
            for (size_t i = 0; i < _syncdims[1]*_syncdims[0]; ++i)
                _sync_r[i] = 0.;
            _syncdims[1]=0;
        }
#else

        _meas  = Matrix<raw>(measdim, 16);
        _meas_r = (float*) &_meas[0];
        _meas_i = _meas_r+1;
        if (syncdim[0]*syncdim[1]>0) {
        	_sync   = Matrix<float>(syncdim, 2);
        	_sync_r = &_sync[0];
        }

#endif
    }

    /**
     * @brief Print parsing results
     */
    void PrintParse () const {
        prtmsg ("     Found %zu lines (data: (Meas: %zu, Noise: %d, Sync: %u))\n",
        		_nlines, _nmeas, 0, _syncdims[1]);
        prtmsg ("       Data dims ( ");
        for (size_t i = 0; i < 15; ++i)
            prtmsg ("%u ", _measdims[i]);
        prtmsg(")\n");
        if (_syncdims[1])
            prtmsg ("       Sync dims ( %ui %ui)\n", _syncdims[0], _syncdims[1] );
        prtmsg ("     Centres: column: %zu, line: %zu, partition: %zu\n", _cent_col, _cent_lin, _cent_par);
        prtmsg ("     Measurement TA: %.2fs, TR: %.2fms\n", _ta, 1000.*_tr);

    }
    
    uint32_t _id;    // Raid file ID
    uint32_t _ndset; // Number of data set in RAID file
    std::vector<EntryHeader> _veh; // Entry headers
    std::vector<uint32_t> _measdims, _syncdims;
    std::vector<size_t> _measdims_sizes; 
    size_t _nmeas, _nlines, _nsync, _cent_par, _cent_lin, _cent_col;
    bool _digested;
    size_t _tend, _tstart;

    float *_meas_r, *_meas_i, *_sync_r, _ta, _tr;
#ifndef USE_IN_MATLAB
        Matrix<raw> _meas, _rtfb;
        Matrix<float> _sync;
#endif
};
}}}}

#endif /* _VD_FILE_HPP_ */
