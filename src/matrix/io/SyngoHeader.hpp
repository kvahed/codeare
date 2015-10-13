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
 *  Created on: Jul, 2014
 *      Author: Kaveh Vahedipour
 */

#ifndef _HEADER_HPP_
#define _HEADER_HPP_

#ifdef _MSC_VER
#  include <cstdint>
#else
#  include <inttypes.h>
#  include <cstdlib>
#endif

#if defined(__GNUC__)
#    define PACKED_MEMBER( __type__, __name__ ) __type__ __name__ __attribute__ ((packed))
#    define PACKED_STRUCT( __type__, __name__ ) __type__ __name__
#else
#    pragma pack(push, pack_header_included)
#    pragma pack(1)
#    define PACKED_MEMBER( __type__, __name__ ) __type__ __name__ 
#    define PACKED_STRUCT( __type__, __name__ ) __type__ __name__
#    pragma pack(pop, pack_header_included)
#endif


#ifndef USE_IN_MATLAB

#include <cassert>
#define prtmsg printf
#define prtwrn printf
#define prterr(x) printf(x); assert(false);
#define evalstr
#define msize_t size_t

#else

#define prtmsg mexPrintf
#define prterr(x) mexErrMsgTxt(x);
#define prtwrn mexWarnMsgTxt
#define evalstr mexEvalString
#define msize_t mwSize

#endif


enum EIM_BIT  {
    ACQEND = 0,
    RTFEEDBACK,
    HPFEEDBACK,
    ONLINE,
    OFFLINE,
    SYNCDATA,
    LASTSCANINCONCAT = 8,
    RAWDATACORRECTION = 10,
    LASTSCANINMEAS,
    SCANSCALEFACTOR,
    SECONDHADAMARPULSE,
    REFPHASESTABSCAN,
    PHASESTABSCAN,
    D3FFT,
    SIGNREV, 
    PHASEFFT,
    SWAPPED, 
    POSTSHAREDLINE,
    PHASCOR,
    PATREFSCAN,
    PATREFANDIMASCAN,
    REFLECT,
    NOISEADJSCAN,
    SHARENOW,
    LASTMEASUREDLINE,
    FIRSTSCANINSLICE,
    LASTSCANINSLICE, 
    TREFFECTIVEBEGIN,
    TREFFECTIVEEND,
};

// Check bit mask
#include <bitset>
static bool bit_set (uint32_t eim, const EIM_BIT eb) {
    std::bitset<sizeof(uint32_t)*8> mask (eim);
    return mask.test(eb);
}

// Check bit mask
#ifdef USE_IN_MATLAB
#include <mex.h>
#endif
static void bit_dump (uint32_t eim) {
    std::bitset<sizeof(uint32_t)*8> mask (eim);
    for (size_t i = 0; i < sizeof(uint32_t)*8; ++i)
        prtmsg("%d", (int)mask.test(i));
    prtmsg("\n");
}

// VB Header
namespace codeare {
namespace matrix  {
namespace io      {
namespace VB {
    typedef struct MeasHeader {
        PACKED_MEMBER(uint32_t, ulDMALength);
        PACKED_MEMBER( int32_t, lMeasUID);
        PACKED_MEMBER(uint32_t, ulScanCounter);
        PACKED_MEMBER(uint32_t, ulTimeStamp);
        PACKED_MEMBER(uint32_t, ulPMUTimeStamp);
        PACKED_MEMBER(uint32_t, aulEvalInfoMask[2]);
        PACKED_MEMBER(uint16_t, ushSamplesInScan);
        PACKED_MEMBER(uint16_t, ushUsedChannels);
        PACKED_MEMBER(uint16_t, sLC[14]);
        PACKED_STRUCT(uint16_t, sCutOff[2]);
        PACKED_MEMBER(uint16_t, ushKSpaceCentreColumn);
        PACKED_MEMBER(uint16_t, ushCoilSelect);
        PACKED_MEMBER(   float, fReadOutOffcentre);
        PACKED_MEMBER(uint32_t, ulTimeSinceLastRF);
        PACKED_MEMBER(uint16_t, ushKSpaceCentreLineNo);
        PACKED_MEMBER(uint16_t, ushKSpaceCentrePartitionNo);
        PACKED_MEMBER(uint16_t, aushIceProgramPara[4]);
        PACKED_MEMBER(uint16_t, aushFreePara[4]);
        PACKED_STRUCT(  float, sSD[7]);
        PACKED_MEMBER(uint16_t, ushChannelId);
        PACKED_MEMBER(uint16_t, ushPTABPosNeg);
    } MeasHeader;
    const size_t MEAS_HEADER_LEN = sizeof(MeasHeader);
}}}}

// VD Header
namespace codeare {
namespace matrix  {
namespace io      {
namespace VD {
    typedef struct MeasHeader {
        PACKED_MEMBER(uint32_t, ulFlagsAndDMALength);
        PACKED_MEMBER( int32_t, lMeasUID);
        PACKED_MEMBER(uint32_t, ulScanCounter);
        PACKED_MEMBER(uint32_t, ulTimeStamp);
        PACKED_MEMBER(uint32_t, ulPMUTimeStamp);
        PACKED_MEMBER(uint16_t, ushSystemType);
        PACKED_MEMBER(uint16_t, ulPTABPosDelay);
        PACKED_MEMBER( int32_t,	lPTABPosX);
        PACKED_MEMBER( int32_t,	lPTABPosY);
        PACKED_MEMBER( int32_t,	lPTABPosZ);
        PACKED_MEMBER(uint32_t,	ulReserved1);
        PACKED_MEMBER(uint32_t, aulEvalInfoMask[2]);
        PACKED_MEMBER(uint16_t, ushSamplesInScan);
        PACKED_MEMBER(uint16_t, ushUsedChannels);
        PACKED_MEMBER(uint16_t, sLC[14]);
        PACKED_STRUCT(uint16_t, sCutOff[2]);
        PACKED_MEMBER(uint16_t, ushKSpaceCentreColumn);
        PACKED_MEMBER(uint16_t, ushCoilSelect);
        PACKED_MEMBER(   float, fReadOutOffcentre);
        PACKED_MEMBER(uint32_t, ulTimeSinceLastRF);
        PACKED_MEMBER(uint16_t, ushKSpaceCentreLineNo);
        PACKED_MEMBER(uint16_t, ushKSpaceCentrePartitionNo);
        PACKED_STRUCT(   float, sSD[7]);
        PACKED_MEMBER(uint16_t, aushIceProgramPara[24]);
        PACKED_MEMBER(uint16_t, aushReservedPara[4]);
        PACKED_MEMBER(uint16_t, ushApplicationCounter);
        PACKED_MEMBER(uint16_t, ushApplicationMask);
        PACKED_MEMBER(uint32_t, ulCRC);
    } MeasHeader;
    const size_t MEAS_HEADER_LEN = sizeof(MeasHeader);

    // Channel header
    typedef struct ChannelHeader {
        PACKED_MEMBER(uint32_t, ulTypeAndChannelLength);
        PACKED_MEMBER( int32_t, lMeasUID);
        PACKED_MEMBER(uint32_t, ulScanCounter);
        PACKED_MEMBER(uint32_t, ulReserved1);
        PACKED_MEMBER(uint32_t, ulSequenceTime);
        PACKED_MEMBER(uint32_t, ulUnused2);
        PACKED_MEMBER(uint16_t, ulChannelId);
        PACKED_MEMBER(uint16_t, ulUnused3);
        PACKED_MEMBER(uint32_t, ulCRC);
    } ChannelHeader;
    const size_t CHANNEL_HEADER_LEN = sizeof(ChannelHeader);
    
    // Entry header
    typedef struct EntryHeader {
        PACKED_MEMBER(uint32_t, MeasID);
        PACKED_MEMBER(uint32_t, FieldID);
        PACKED_MEMBER(uint64_t, MeasOffset);
        PACKED_MEMBER(uint64_t, MeasLen);
        char PatientName [64];
        char ProtocolName[64];
    } EntryHeader;
    const size_t ENTRY_HEADER_LEN = sizeof(EntryHeader);

}}}}

// Sync header length
const size_t SYNC_HEADER_SIZE = 64;

#endif   /* _HEADER_HPP_ */
