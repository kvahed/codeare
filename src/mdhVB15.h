/**
 * VB15A MDH
 */ 

#ifndef __MDH_H__
#define __MDH_H__

typedef struct {
	
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


// Slice vectors
typedef struct {
	
    float  flSag;
    float  flCor;
    float  flTra;

} sVector;

// Slice data
typedef struct {
	
	sVector  sSlicePosVec;                    /* slice position vector        */
	float    aflQuaternion[4];                /* rotation matrix as quaternion*/
	
} sSliceData;                               /* sizeof : 28 byte             */

// Cut off
typedef struct {

	unsigned short  ushPre;               /* write ushPre zeros at line start */
	unsigned short  ushPost;              /* write ushPost zeros at line end  */

} sCutOffData;

// MDH
typedef struct {

	unsigned int   ulDMALength;                // DMA length [bytes] must be                     4 
	int            lMeasUID;                   // measurement user ID                            4     
	unsigned int   ulScanCounter;              // scan counter [1...]                            4
	unsigned int   ulTimeStamp;                // time stamp [2.5 ms ticks since 00:00]          4
	unsigned int   ulPMUTimeStamp;             // PMU time stamp [2.5 ms since last trigge]      4
	unsigned int   aulEvalInfoMask[2];         // evaluation info mask field                     8
	unsigned short ushSamplesInScan;           // # of samples acquired in scan                  2
	unsigned short ushUsedChannels;            // # of channels used in scan                     2 =  32
	sLoopCounter   sLC;                        // loop counters                                 28 =  60
	sCutOffData    sCutOff;                    // cut-off values                                 4           
	unsigned short ushKSpaceCentreColumn;      // centre of echo                                 2
	unsigned short ushCoilSelect;              // for swapping                                   2
	float          fReadOutOffcentre;          // ReadOut offcenter value                        4
	unsigned int   ulTimeSinceLastRF;          // Sequence time stamp since last RF pulse        4
	unsigned short ushKSpaceCentreLineNo;      // number of K-space centre line                  2
	unsigned short ushKSpaceCentrePartitionNo; // number of K-space centre partition             2
	unsigned short aushIceProgramPara[4];      // free parameter for IceProgram                  8 =  88
	unsigned short aushFreePara[4];            // free parameter                         4 * 2 = 8   
	sSliceData     sSD;                        // Slice Data                                    28 = 124
	unsigned short ushChannelId;               // channel Id must be the last parameter          2
	unsigned short ushPTABPosNeg;              // negative, absolute PTAB position in [0.1 mm]   2

} sMDH;                                        // total length:                                  0 = 128 bytes

#endif   /* MDH_H */
