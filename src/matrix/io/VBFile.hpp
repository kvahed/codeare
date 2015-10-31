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

#ifndef _VB_FILE_HPP_
#define _VB_FILE_HPP_

#include "SyngoFile.hpp"
#include "Workspace.hpp"
#ifdef USE_IN_MATLAB
  #include "waitmex.h"
#endif

#include "SimpleTimer.hpp"

namespace codeare {
namespace matrix  {
namespace io      {

namespace VB {
    class VBFile : public SyngoFile {
        
    public:
#ifdef USE_MATLAB
        VBFile (const std::string fname, int nlhs = 0, mxArray *lhs[] = 0) :
            SyngoFile(fname, nlhs, lhs), _nnoise (0), _nmeas(0), _nctnorm(0), _data_len(0),
            _meas_r(0), _meas_i(0), _noise_r(0), _noise_i(0), _sync_r(0), _lines(0),
            _nrtfb(0), _rtfb_r(0), _rtfb_i(0), _meas(0) {
            _file.seekg (0, _file.end);
            _data_len = _file.tellg();
            _file.seekg (0, _file.beg);
            _file.read ((char*)&_header_len, sizeof(uint32_t));
            _file.seekg(0);
            _measdims = std::vector<uint32_t>(16);
            _rtfbdims = std::vector<uint32_t>(16);
            _syncdims = std::vector<uint32_t>( 2);
            prtmsg("     VB file contains %.1fGB of raw data starting at position %d.\n", 1e-9f*iGB*_data_len, _header_len);
        }
#else
        VBFile (const std::string& fname, const IOMode mode, const Params& params, const bool verbosity) :
        	SyngoFile(fname, mode, params, verbosity), _nnoise (0), _nmeas(0), _nctnorm(0), _data_len(0),
			_meas_r(0), _meas_i(0), _noise_r(0), _noise_i(0), _sync_r(0), _lines(0), _nrtfb(0),
			_rtfb_r(0), _rtfb_i(0), _digested(false), _tend(0), _tstart(0) {
            _file.seekg (0, _file.end);
            _data_len = _file.tellg();
            _file.seekg (0, _file.beg);
            _file.read ((char*)&_header_len, sizeof(uint32_t));
            _file.seekg(0);
            _measdims = std::vector<uint32_t>(16);
            _rtfbdims = std::vector<uint32_t>(16);
            _syncdims = std::vector<uint32_t>( 2);
            prtmsg("     VB file contains %.1fGB of raw data starting at position %d.\n", 1e-9f*iGB*_data_len, _header_len);
            Digest();
            wspace.Add<cxfl>("meas", _meas);
            wspace.Add<cxfl>("rtfb", _rtfb);
            wspace.Add<float>("sync", _sync);
        }
#endif

        virtual ~VBFile () {
            _meas_r = 0;
            _meas_i = 0;
            _sync_r = 0;
            _syncdims.clear();
            _measdims.clear();
            _noisedims.clear();
        }
        
        
        virtual void Digest () {

            MeasHeader mh;
            size_t i = 0;
            char* wbm = const_cast<char*>("reading ...");
            SimpleTimer st;
            
            _file.seekg (_header_len);
            
            if (!_meas_r)
                prtmsg ("   Parsing ... \n");
            else {
                prtmsg ("   Reading ... \n");
                //h = waitbar_create (0, wbm) ;
            }
            
            while (true) {

                _file.read ((char*)&mh, MEAS_HEADER_LEN);
                fflush(stdout);

                if (bit_set(mh.aulEvalInfoMask[0], ACQEND)) {
                    prtmsg ("     Hit ACQEND\n");
                    break;
                } else if (bit_set(mh.aulEvalInfoMask[0], SYNCDATA)) {
                    size_t tmp1 = _file.tellg(), tmp2;
                    ParseSync (mh);
                    tmp2 = _file.tellg();
                    tmp1 = ((tmp2-tmp1)%32);
                    if (tmp1)
                        _file.seekg (tmp2 + 32-tmp1);
                } else if (bit_set(mh.aulEvalInfoMask[0], NOISEADJSCAN)) {
                    ParseNoise (mh);
                } else if (bit_set(mh.aulEvalInfoMask[0], RTFEEDBACK)) {
                    if (_nrtfb == 0 && !_rtfb_r) {
                        _rtfbdims[0] = mh.ushSamplesInScan;
                        _rtfbdims[1] = mh.ushUsedChannels;
                    }
                    ParseFeedback (mh);
                } else if (bit_set(mh.aulEvalInfoMask[1], ONLINE)) {
                    ParseCTNormalize (mh);
                } else {
                    if (_nmeas == 0 && !_meas_r) {
                        _measdims[0] = mh.ushSamplesInScan;
                        _measdims[1] = mh.ushUsedChannels;
                    }
                    ParseMeas (mh);
                }
                if (_meas_r && i%1000==0) {
                    //waitbar_update ((double)i/(double)_lines, h, wbm);
                }
                
                ++i;
            }
            _lines = i;
            if (_meas_r) {
                //waitbar_destroy(h);
            }
            prtmsg ("     done - wtime %s", st.Format().c_str());
            
            if (!_meas_r) {
                _measdims = _raise_one (_measdims);
                _rtfbdims = _raise_one (_rtfbdims);
                _measdims[5] = (_measdims[5]-_cent_par)*2;
                _ta = 2.5e-3*(_tend-_tstart);
                wspace.PSet("TA", _ta);
                _tr = _ta/_nmeas;
                prtmsg ("     Found %d lines (data: (Meas: %d, Noise: %d, Sync: %d))\n", i, _nmeas, _nnoise, _syncdims[1]);
                prtmsg ("     Centres: column: %zu, line: %zu, partition: %zu\n", _cent_col, _cent_lin, _cent_par);
                prtmsg ("     Measurement TA: %.2fs, TR: %.2fms\n", _ta, 1000.*_tr);
                prtmsg ("       Data dims ( ");
                for (size_t i = 0; i < 15; ++i)
                    prtmsg ("%d ", _measdims[i]);
                prtmsg(")\n");
                prtmsg ("       RTFB dims ( ");
                for (size_t i = 0; i < 15; ++i)
                    prtmsg ("%d ", _rtfbdims[i]);
                prtmsg(")\n");
                prtmsg ("       Sync dims ( %d %d)\n", _syncdims[0], _syncdims[1] );
                evalstr("drawnow;");
                this->Allocate();
                this->Digest();
            }
            _digested=true;
            
        }
        
        template<class T> Matrix<T> Read (const std::string& data_name) {

        	if (!_digested) {
        		this->Digest();
                wspace.Add<cxfl>("meas", _meas);
                wspace.Add<cxfl>("rtfb", _rtfb);
                wspace.Add<float>("sync", _sync);
            }

        	return Matrix<T>();
        }


    private:
        
        void Allocate () {
            msize_t measdim[16] = {
                _measdims[ 0], _measdims[ 1], _measdims[ 2], _measdims[ 3], _measdims[ 4], _measdims[ 5], _measdims[ 6], _measdims[ 7],
                _measdims[ 8], _measdims[ 9], _measdims[10], _measdims[11], _measdims[12], _measdims[13], _measdims[14], _measdims[15]};
            msize_t rtfbdim[16] = {
                _rtfbdims[ 0], _rtfbdims[ 1], _rtfbdims[ 2], _rtfbdims[ 3], _rtfbdims[ 4], _rtfbdims[ 5], _rtfbdims[ 6], _rtfbdims[ 7],
                _measdims[ 8], _rtfbdims[ 9], _rtfbdims[10], _rtfbdims[11], _rtfbdims[12], _rtfbdims[13], _rtfbdims[14], _rtfbdims[15]};
            msize_t syncdim[2] = {_syncdims[0], _syncdims[1]};
#ifdef USE_IN_MATLAB
            // Measurement data
            _lhs[0] = mxCreateNumericArray (16, measdim, mxSINGLE_CLASS, mxCOMPLEX);
            _meas_r = (float*) mxGetData(_lhs[0]);
            _meas_i = (float*) mxGetImagData(_lhs[0]);
            _nmeas = 0;
            
            // Realtime
            if (_nlhs >= 2) {
                _lhs[1] = mxCreateNumericArray (16, rtfbdim, mxSINGLE_CLASS, mxCOMPLEX);
                _rtfb_r = (float*) mxGetData(_lhs[1]);
                _rtfb_i = (float*) mxGetImagData(_lhs[1]);
            }
            
            // Sync data
            if (_nlhs >= 3) {
                _lhs[2] = mxCreateNumericArray (2, syncdim, mxSINGLE_CLASS, mxREAL);
                _sync_r = (float*) mxGetData(_lhs[2]);
            }
#else
            _meas  = Matrix<raw>(measdim, 16);
            _meas_r = (float*) &_meas[0];
            _meas_i = _meas_r+1;
            if (rtfbdim[0]*rtfbdim[1]>0) {
                _rtfb   = Matrix<raw>(rtfbdim, 16);
                _rtfb_r = (float*) &_rtfb[0];
                _rtfb_i = _rtfb_r+1;
            }
            if (syncdim[0]*syncdim[1]>0) {
            	_sync   = Matrix<float>(syncdim, 2);
            	_sync_r = &_sync[0];
            }
#endif
            _nmeas = 0;
            _nrtfb = 0;
            _syncdims[1]=0;
        }
        
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
        
        inline void ParseNoise (const MeasHeader& mh) {
            std::vector<std::complex<float> > noise (mh.ushSamplesInScan);
            _file.read((char*)&noise[0], mh.ushSamplesInScan*sizeof(std::complex<float>));
            if (!_noise_r)
                _nnoise++;
        }
        
        inline void ParseCTNormalize (const MeasHeader& mh) {
            std::vector<std::complex<float> > buf (mh.ushSamplesInScan);
            _file.read((char*)&buf[0], mh.ushSamplesInScan*sizeof(std::complex<float>));
        }
        
        inline void ParseMeas (const MeasHeader& mh) {
            size_t pos (_file.tellg());
            if (!_meas_r) {
            	if (_nmeas == 0 && mh.ushChannelId == 0) {
            		_tstart   = mh.ulTimeStamp;
            		_cent_par = mh.ushKSpaceCentrePartitionNo;
            		_cent_lin = mh.ushKSpaceCentreLineNo;
            		_cent_col = mh.ushKSpaceCentreColumn;
            	} else if (bit_set(mh.aulEvalInfoMask[0], LASTSCANINMEAS) && mh.ushChannelId == 0) {
            		_tend = mh.ulTimeStamp;
            	}
                _measdims = _max (_measdims, mh.sLC);
                _file.seekg (pos + _measdims[0]*sizeof(std::complex<float>));
            } else {
#ifndef USE_IN_MATLAB
               	_file.read ((char*)&_meas(0, mh.ushChannelId, mh.sLC[0], mh.sLC[1], mh.sLC[2],
               			mh.sLC[3], mh.sLC[4], mh.sLC[5], mh.sLC[6], mh.sLC[7], mh.sLC[8], mh.sLC[9],
						mh.sLC[10], mh.sLC[11], mh.sLC[12], mh.sLC[13]),
						mh.ushSamplesInScan*sizeof(std::complex<float>));
#else
                std::vector<std::complex<float> >buf (mh.ushSamplesInScan);
                _file.read((char*)&buf[0], mh.ushSamplesInScan*sizeof(std::complex<float>));
                for (size_t i = 0; i < _measdims[0]; ++i) {
                    _meas_r[i+_nmeas*_measdims[0]] = std::real(buf[i]);
                    _meas_i[i+_nmeas*_measdims[0]] = std::imag(buf[i]);
                }
#endif
            }
            _nmeas++;
        }
        
        inline void ParseFeedback (const MeasHeader& mh) {
            size_t pos (_file.tellg());
            if (!_rtfb_r) {
                _rtfbdims = _max (_rtfbdims, mh.sLC);
                _file.seekg (pos + _rtfbdims[0]*sizeof(std::complex<float>));
            } else {
#ifndef USE_IN_MATLAB
                _file.read((char*)&_rtfb[_nrtfb*_rtfbdims[0]], mh.ushSamplesInScan*sizeof(std::complex<float>));
#else
                std::vector<std::complex<float> >buf (mh.ushSamplesInScan);
                _file.read((char*)&buf[0], mh.ushSamplesInScan*sizeof(std::complex<float>));
                for (size_t i = 0; i < _rtfbdims[0]; ++i) {
                    _rtfb_r[i+_nrtfb*_rtfbdims[0]] = std::real(buf[i]);
                    _rtfb_i[i+_nrtfb*_rtfbdims[0]] = std::imag(buf[i]);
                }
#endif
            }
            _nrtfb++;
        }

        bool _digested;
        size_t _data_len, _cent_par, _cent_lin, _cent_col;

        unsigned _nmeas, _nrtfb, _nnoise, _lines, _nctnorm;
        float *_meas_r, *_meas_i, *_rtfb_r, *_rtfb_i, *_noise_r, *_noise_i, *_sync_r, _ta, _tr;
        size_t _tend, _tstart;
        std::vector<uint32_t> _syncdims, _measdims, _rtfbdims, _noisedims;
#ifndef USE_IN_MATLAB
        Matrix<raw> _meas, _rtfb, _ctnorm;
        Matrix<float> _sync;
#else
        void _meas;
#endif
        
    };
}
}}}
#endif /* _VD_FILE_HPP_ */
