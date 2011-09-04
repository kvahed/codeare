/*
 *mdhVB15.h 
 *
 * header file containing the definition of the structs of the MDHs stored in VB15 .dat files
 *THIS IS NOT THE ORIGINAL MDH.H FILE USED BY SIEMENS. THIS IS A MODIFIED VERSION (S.Romanzetti)
 * 
 */
/*--------------------------------------------------------------------------*/
/* Include control                                                          */
/*--------------------------------------------------------------------------*/
#ifndef MDH_H
#define MDH_H


typedef struct
{
  unsigned short  ushLine;                  /* line index                   */
  unsigned short  ushAcquisition;           /* acquisition index            */
  unsigned short  ushSlice;                 /* slice index                  */
  unsigned short  ushPartition;             /* partition index              */
  unsigned short  ushEcho;                  /* echo index                   */	
  unsigned short  ushPhase;                 /* phase index                  */
  unsigned short  ushRepetition;            /* measurement repeat index     */
  unsigned short  ushSet;                   /* set index                    */
  unsigned short  ushSeg;                   /* segment index  (for TSE)     */
  unsigned short  ushIda;                   /* IceDimension a index         */
  unsigned short  ushIdb;                   /* IceDimension b index         */
  unsigned short  ushIdc;                   /* IceDimension c index         */
  unsigned short  ushIdd;                   /* IceDimension d index         */
  unsigned short  ushIde;                   /* IceDimension e index         */
} sLoopCounter;   

/*--------------------------------------------------------------------------*/
/*  Definition of slice vectors                                             */
/*--------------------------------------------------------------------------*/

typedef struct
{
    float  flSag;
    float  flCor;
    float  flTra;
} sVector;

typedef struct
{
  sVector  sSlicePosVec;                    /* slice position vector        */
  float    aflQuaternion[4];                /* rotation matrix as quaternion*/
} sSliceData;                               /* sizeof : 28 byte             */


/*--------------------------------------------------------------------------*/
/*  Definition of cut-off data                                              */
/*--------------------------------------------------------------------------*/
typedef struct
{
  unsigned short  ushPre;               /* write ushPre zeros at line start */
  unsigned short  ushPost;              /* write ushPost zeros at line end  */
} sCutOffData;



/*--------------------------------------------------------------------------*/
/*  Definition of measurement data header                                   */
/*--------------------------------------------------------------------------*/
typedef struct
{
  unsigned int  ulDMALength;                  // DMA length [bytes] must be                         4 byte                                               // first parameter                        
  int            lMeasUID;                     // measurement user ID                                4     
  unsigned int   ulScanCounter;                // scan counter [1...]                                4
  unsigned int   ulTimeStamp;                  // time stamp [2.5 ms ticks since 00:00]              4
  unsigned int   ulPMUTimeStamp;               // PMU time stamp [2.5 ms ticks since last trigger]   4
  unsigned int   aulEvalInfoMask[2];           // evaluation info mask field                         8
  unsigned short ushSamplesInScan;             // # of samples acquired in scan                     2
  unsigned short ushUsedChannels;              // # of channels used in scan                        2   =32
  sLoopCounter   sLC;                          // loop counters                                    28   =60
  sCutOffData    sCutOff;                      // cut-off values                                    4           
  unsigned short ushKSpaceCentreColumn;        // centre of echo                                    2
  unsigned short ushCoilSelect;                     // for swapping                                      2
  float          fReadOutOffcentre;            // ReadOut offcenter value                           4
  unsigned int   ulTimeSinceLastRF;            // Sequence time stamp since last RF pulse           4
  unsigned short ushKSpaceCentreLineNo;        // number of K-space centre line                     2
  unsigned short ushKSpaceCentrePartitionNo;   // number of K-space centre partition                2
  unsigned short aushIceProgramPara[4];        // free parameter for IceProgram               8   =88
  unsigned short aushFreePara[4];              // free parameter                          4 * 2 =   8   
  sSliceData     sSD;                          // Slice Data                                 28  =124
  unsigned short ushChannelId;                  // channel Id must be the last parameter             2
  unsigned short ushPTABPosNeg;                  // negative, absolute PTAB position in [0.1 mm]     2
                                                                   // (automatically set by PCI_TX firmware)
} sMDH;                                        // total length: 32 * 32 Bit (128 Byte)            128


#endif   /* MDH_H */
