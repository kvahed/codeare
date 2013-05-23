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

#ifndef __IO_HPP__
#define __IO_HPP__

enum io_strategy {

	HDF5,
	MATLAB,
	NIFTI,
	SYNGO,
	CDF,
	PRIMITIVE

};

#include "Matrix.hpp"
#include "Toolbox.hpp"
#include "Tokenizer.hpp"
#include "Waitbar.hpp"
#include "tinyxml/tinyxml.h"
#include "tinyxml/xpath_static.h"
#include "mdhVB15.h"
#include "Algos.hpp"

#ifdef HAVE_NIFTI1_IO_H
#include <nifti1_io.h>
#endif

#ifdef HAVE_H5CPP_H
#include <H5Cpp.h>
using namespace H5;
#endif

#ifdef HAVE_CDF_H
#include <cdf.h>
#endif


#include <limits>
#include <map>
#include <sstream>
#include <iomanip>
#ifdef HAVE_CXXABI_H
    #include <cxxabi.h>
#endif

using namespace TinyXPath;
using namespace std;


enum dtype {

	RLFL,
	RLDB,
	CXFL,
	CXDB,
	LONG,
	SHRT

};


/**
 * @brief           Verbose wrapper around fwrite
 * 
 * @param  d        Data repository to write from
 * @param  sz       Size of individual elements
 * @param  n        Number of elements
 * @param  f        File handle
 * @param  desc     Description
 * @return          Success
 */
inline static bool 
mwrite (const void* d, const size_t sz, const size_t n, FILE* f, std::string desc) {

	if (size_t l = fwrite (d, sz, n, f) != n) {

		printf("File write error - %s: %li != %li!\n", desc.c_str(), l, n);
		fclose (f);
		return false;
		
	}

	return true;

} 
	

/**
 * @brief            Primitive dump of data.<br> The data id dumped column major into a file.
 *
 * @param   M        Matrix
 * @param   fname    File name       
 * @return           Success
 */
template <class T> inline static bool
PRDump (const Matrix<T>& M, const string fname) {
	
	FILE *f;
	size_t n = numel(M);
	int dt;
	size_t ns, l;

	if      (typeid(T) == typeid(float))  dt = RLFL;
	else if (typeid(T) == typeid(double)) dt = RLDB;
	else if (typeid(T) == typeid(cxfl))   dt = CXFL;
	else if (typeid(T) == typeid(cxdb))   dt = CXDB;
	else if (typeid(T) == typeid(long))   dt = LONG;
	else if (typeid(T) == typeid(short))  dt = SHRT;
	
	if ((f = fopen(fname.c_str(), "wb"))==NULL) {
		printf("Cannot open %s file for writing.\n", fname.c_str());
		return false;
	}

	// Dump type
	if (!mwrite(&dt, sizeof(int), 1, f, "data type"))
		return false;
	
	// Dump dimensions
	if (!mwrite((const void*) M.Dim(), sizeof(size_t), INVALID_DIM, f, "dimensions"))
		return false;
	
	// Dump resolutions
	if (!mwrite((const void*) M.Res(), sizeof(float), INVALID_DIM, f, "resolutions"))
		return false;
	
	// Size of name and name
	ns = strlen(M.GetClassName());
	if (!mwrite(&ns, sizeof(size_t), 1, f, "name length"))
		return false;
	
	// Dump name
	if (!mwrite((const void*) M.GetClassName(), sizeof(char), ns, f, "name"))
		return false;
	
	// Dump data
	if (!mwrite((const void*) M.Memory(), sizeof(T), n, f, "data"))
		return false;
	
	fclose(f);
	
	return true;
	
}


/**
 * @brief           Verbose wrapper around fread
 * 
 * @param  d        Place holder to read into
 * @param  sz       Size of individual elements
 * @param  n        Number of elements
 * @param  f        File handle
 * @param  desc     Description
 * @return          Success
 */
inline static bool 
mread (void* d, const size_t sz, const size_t n, FILE* f, const std::string desc) {

	if (size_t l = fread (d, sz, n, f) != n) {

		printf("File read error - %s: %li != %li!\n", desc.c_str(), l, n);
		fclose (f);
		return false;
		
	}

	return true;

} 
	

/**
 * @brief            Primitive dump of data.<br> The data id dumped column major into a file.
 *
 * @param   M        Matrix
 * @param   fname    File name       
 * @return           Success
 */
template <class T> inline static bool
PRRead (Matrix<T>& M, const string fname) {
	
	FILE *f;
	int dt;
	size_t  ns, n;
	std::vector<size_t> dims (INVALID_DIM,1);
	float res[INVALID_DIM];
	char* name;
	
	if ((f = fopen(fname.c_str(), "rb"))==NULL) {
		printf("Cannot open %s file for reading.\n", fname.c_str());
		return false;
	}

	// Read type
	if (!mread ( &dt, sizeof(   int),           1, f, "data type")) return false;

	// Matrix and data type must fit as of now.
	if ((typeid(T) == typeid(cxfl) && dt == RLFL) ||
		(typeid(T) == typeid(cxfl) && dt == RLDB) ||
		(typeid(T) == typeid(cxfl) && dt == CXFL) ||
		(typeid(T) == typeid(cxfl) && dt == CXDB) ||
		(typeid(T) == typeid(cxfl) && dt == LONG) ||
		(typeid(T) == typeid(cxfl) && dt == SHRT)) {
	
		// Read dimensions and allocate matrix
		if (!mread (&dims[0], sizeof(size_t), INVALID_DIM, f, "dimensions")) return false;
		M = Matrix<T>(dims);
		n = numel(M);
		
		// Read resolutions and assign
		if (!mread ( res, sizeof( float), INVALID_DIM, f, "resolutions")) return false;
		for (size_t i = 0; i < INVALID_DIM; i++)
			M.Res(i) = res[i];
		
		// Name 
		if (!mread ( &ns, sizeof(size_t),           1, f, "name length")) return false;
		name = new char [ns];
		if (!mread (name,   sizeof(char),          ns, f, "name")) return false;
		M.SetClassName(name);
		
		// Read data
		if (!mread (&M[0],     sizeof(T),           n, f, "data")) return false;

		//Close and clean up;
		fclose(f);
		delete name;

	} else {

		fclose (f);
		return false;

	}
		
	return true;
	
}
	

/**
 * @brief         Prepare matrix for Syngo MR read.
 *
 * @param  M      Matrix
 * @param  fname  File name
 * @return        Success
 */
template <class T> static bool
RSAdjust (Matrix<T>& M, const std::string& fname) {
	
	size_t dimk;
	size_t dimv;
	
	size_t ndims [INVALID_DIM];
	float  nress [INVALID_DIM];
	
	// XML file name, run cmd and repository
	stringstream    xmlf;
	stringstream    cmd;
	
	// XML objects
	TiXmlDocument*       doc = new TiXmlDocument();
	TiXmlElement*       meta;
	TiXmlNode*           vol;
	
	// Dimensions in XProtocol
	
	map < string, size_t > dims;
	dims["NImageCols"] =  0; dims["NLinMeas"] =  1; dims["NSlcMeas"] =  2; dims["NParMeas"] =  3;
	dims[  "NEcoMeas"] =  4; dims["NPhsMeas"] =  5; dims["NRepMeas"] =  6; dims["NSetMeas"] =  7;
	dims[  "NSegMeas"] =  8; dims[  "RawCha"] =  9; dims["NIdaMeas"] = 10; dims["NIdbMeas"] = 11;
	dims[  "NIdcMeas"] = 12; dims["NIddMeas"] = 13; dims["NIdeMeas"] = 14; dims["NAveMeas"] = 15; 
	
	string channels  = "RawCha";
	string doublevals [2] = {"RoFOV", "PeFOV"};
	string slcthickn = "SliceThickness";
	
	// Create XML output from XProt
	cmd << "/usr/local/bin/convert.pl ";
	cmd << fname;
	size_t ec = system (cmd.str().c_str());
	
	// Output filename
	xmlf << fname;
	xmlf << ".xml";
	
	// Read XML and go to rawobjectprovider for dimensions
	doc->LoadFile (xmlf.str().c_str());
	meta = doc->RootElement();
	
	TiXmlNode* rootfunctor    = TinyXPath::XNp_xpath_node 
		(meta, "/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']");
	
	// Dimensions
	vol  = rootfunctor->FirstChild ("ParamLong");
	map < string, size_t >::iterator it;
	string key;
	
	do {
		
		key = ((TiXmlElement*)vol)->Attribute("name");
		it = dims.find (key);
		
		if (it != dims.end()) {
			
			dimk = dims[key];
			dimv = atoi( ((TiXmlElement*)vol)->Attribute("value"));
			
			if (dimk == 0) dimv *= 2;
			if (dimv == 0) dimv = 1;
			
			ndims[dimk] = dimv;
			
		}
		
	} while ((vol = rootfunctor->IterateChildren ("ParamLong", vol))!=NULL);
	
	nress[0] = (float)TinyXPath::d_xpath_double
		(meta, "/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamDouble[@name='RoFOV']/@value") / ndims[0] * 2;
	nress[1] = (float)TinyXPath::d_xpath_double 
		(meta, "/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamDouble[@name='PeFOV']/@value") / ndims[1];
	nress[3] = (float)TinyXPath::d_xpath_double
		(meta, "/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online1']/ParamFunctor[@name='root']/ParamArray[@name='SliceThickness']/ParamDouble[1]/Precision[1]/@value") / ndims[3];
	
	ndims[9] = i_xpath_int 
		(meta, "/Meta/XProtocol[1]/ParamMap[1]/ParamMap[1]/Pipe[1]/PipeService[@name='Online2']/ParamFunctor[@name='rawobjprovider']/ParamLong[@name='RawCha']/@value"); 
	
	M = Matrix<T> (ndims, nress);
	
	delete doc;
	return true;
	
}


/**
 * @brief          Read from Syngo MR .dat file
 * 
 * @param  M       Matrix to hold the data
 * @param  fname   File name
 * @param  version Syngo MR version 
 * @return         Success
 */
template <class T> static bool 
RAWRead (Matrix<T>& M, const std::string& fname, const std::string& version) {
	
	// Get size information from XProtocol and resize 
	ticks  tic = getticks();
	
	RSAdjust(M, fname);
	
	printf ("%s: %s\n", fname.c_str(), DimsToCString(M));
	
	FILE*         f;
	sMDH*         mdh;
	unsigned      l      = 0;
	size_t        nscans = (numel(M) / size(M,0));
	unsigned      start;
	unsigned      size, read;
	
	// Assess data size
	f = fopen (fname.c_str(), "rb");
	read = fread (&l, sizeof(unsigned), 1, f);
	fseek (f,     l, SEEK_SET);
	start = ftell(f);
	fseek (f,    -1, SEEK_END);
	size  = ftell(f);
	fseek (f, start, SEEK_SET);
	
	//printf ("Found %.1f MB of header and data.\n", (float)size / MB);
	stringstream pre;
	
	float bsize = (float) size / (1024.0*1024.0);
	pre << "  Reading " << bsize << "MB ... ";
	
	// Headers
	mdh  = (sMDH*) malloc (nscans * sizeof(sMDH));
	size_t n    = 0.0;
	char   mask[8];
	bool   done = false;
	
	const std::string post = "";
	
	size_t i   = 0;
	
	for (i = 0; i < nscans; i++) {
		
		if (i % 250 == 0 || i == nscans-1)
			ProgressBar (pre.str(), post, (short)round((float)(i+1)/(float)nscans*100.0));
		
		read = fread (&mdh[i], sizeof(sMDH), 1, f);
		
		long m = 0;
		memcpy (&m, mdh[i].aulEvalInfoMask, 2 * sizeof (unsigned));
		if (m == 0x00000001) {
			ProgressBar (pre.str(), post, 100);
			done = true;
				break;
		}
		
		n = mdh[i].sLC.ushLine       * M.Dim(0) +
			mdh[i].sLC.ushSlice      * M.Dim(0) * M.Dim(1) +
			mdh[i].sLC.ushPartition  * M.Dim(0) * M.Dim(1) * M.Dim(2) +
			mdh[i].sLC.ushEcho       * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) +
			mdh[i].sLC.ushPhase      * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) * M.Dim(4) +
			mdh[i].sLC.ushRepetition * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) * M.Dim(4) * M.Dim(5) +
			mdh[i].sLC.ushSet        * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) * M.Dim(4) * M.Dim(5) * M.Dim(6) +
			mdh[i].sLC.ushSeg        * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) * M.Dim(4) * M.Dim(5) * M.Dim(6) * M.Dim(7) +
			mdh[i].ushChannelId      * M.Dim(0) * M.Dim(1) * M.Dim(2) * M.Dim(3) * M.Dim(4) * M.Dim(5) * M.Dim(6) * M.Dim(7) * M.Dim(8);
		
		read = fread (&M[n], sizeof (complex<float>), M.Dim(0), f);
		
	}
	
	fclose (f);
	
	float e = ((float)elapsed(getticks(),tic) / (float)Toolbox::Instance()->ClockRate());
	printf ("\n  done. Loaded %i records (%f MB/s).\n", (int)i, bsize/e);
	
	return true;
	
}



/**
 * @brief          Dump matrix to NIFTI file
 * 
 * @param  M       Matrix
 * @param  fname   File name 
 * @return         Success
 */
template <class T> static bool
NIDump (const Matrix<T>& M, const string fname) {
	
	if (fname != "") {
		
#ifdef HAVE_NIFTI1_IO_H
		
		Matrix<T>      tmp = M;
		tmp = squeeze(tmp);
		
		size_t            l   = fname.length();
		
		nifti_1_header header;
		header.sizeof_hdr = 348;
		header.dim[0] = ndims(tmp);
		header.pixdim[0] = ndims(tmp);
		
		if (ndims(tmp) > 7) {
			printf ("Cannot dump more than 8 dimensions to NIFTI FILE\n.");
			return false;
		}
		
		for (size_t i = 0; i < 7; i++) {
			header.dim   [i+1] = tmp.Dim(i);
			header.pixdim[i+1] = tmp.Res(i);
		}
		
		if      (typeid(T) == typeid(cxfl))
			header.datatype = 32;
		else if (typeid(T) == typeid(double))
			header.datatype = 64;
		else if (typeid(T) == typeid(short))
			header.datatype = 256;
			
		nifti_image* ni = nifti_convert_nhdr2nim(header, NULL);
		
		ni->nifti_type = 1;
		
		// Single nii.gz file
		ni->fname = (char*) calloc(1,l); 
		strcpy(ni->fname,fname.c_str());
		ni->iname = (char*) calloc(1,l); 
		strcpy(ni->iname,fname.c_str());
		
		ni->data = (void*) malloc (numel(M) * sizeof (T));
		memcpy (ni->data, M.Memory(), numel(M) * sizeof (T));
		
		nifti_image_write (ni);
		nifti_image_free (ni); 
		
		return true;
		
#else 
		
		printf ("NIFTI IO ERROR - Didn't dump nothin'");
		return false;
		
#endif
		
	} else {
		
		printf ("NIFTI IO ERROR - Empty file name");
		return false;
		
	}
	
}
	

/**
 * @brief          Read matrix from NIFTI file
 * 
 * @param  M       Matrix
 * @param  fname   File name 
 * @return         Success
 */
template <class T> static bool
NIRead (Matrix<T>& M, const string fname) {
	
	std::vector<size_t> ndims (INVALID_DIM,1);
	std::vector<float>  nress (INVALID_DIM,1.0);
	size_t i = 0;

	if (fname != "") {
		
#ifdef HAVE_NIFTI1_IO_H
		
		nifti_image* ni = nifti_image_read (fname.c_str(), 1);
		
		if (ni == NULL) 
			return false;
		
		for (; i < ni->dim[0]; i++)
			if (ni->dim[i+1] > 1) {
				ndims[i] = ni->dim[i+1];
				nress[i] = ni->pixdim[i+1];
			}
		
		for (; i < INVALID_DIM; i++) {
			ndims[i] = 1;
			nress[i] = 1.0;
		}
		
		M = Matrix<T>(ndims, nress);
		
		if ((ni->datatype == 16 || ni->datatype == 64) && typeid(T) == typeid(double)) {
			if (ni->datatype == 64)
				memcpy (&M[0], ni->data,numel(M)*sizeof(T));
			else 
				for (size_t i = 0; i <numel(M); i++ )
					M[i] = ((float*)ni->data)[i];
		} else if ((ni->datatype == 32 || ni->datatype == 1792) && typeid(T) == typeid(cxfl)) {
			if (ni->datatype == 32)
				memcpy (&M[0], ni->data,numel(M)*sizeof(T));
			else 
				for (size_t i = 0; i <numel(M); i++) {
					float f[2] = {((double*)ni->data)[2*i], ((double*)ni->data)[2*i+1]};
					memcpy(&M[i], f, 2 * sizeof(float));
				}
		} else if ((ni->datatype == 256 || ni->datatype == 4) && typeid(T) == typeid(short)) {
			if (ni->datatype == 256 || ni->datatype == 4)
				memcpy (&M[0], ni->data, SizeInRAM(M));
			
		} else {
			printf (" Unsupported data type %i!", ni->datatype);
			return false;
		}
		
		nifti_image_free (ni);
			
		return true;
		
#else 

		printf ("NIFTI IO ERROR - Missing NIFTI libraries/headers");
		return false;
			
#endif
			
	} else {
		
		printf ("NIFTI IO ERROR - Empty file name");
		return false;
		
	}
	
}
	
	
/**
 * @brief          Read matrix to CDF file.<br/> !!! FUNCTION NOT READY YET !!!
 * 
 * @param  M       Matrix
 * @param  fname   File name 
 * @param  dname   Data name
 * @param  dloc    Data location
 * @return         Success
 */
template <class T> static bool 
CDFDump (const Matrix<T>& M, const string fname, const string dname, const string dloc) {
	
#ifdef HAVE_CDF_H
		
	CDFid     id;                // CDF identifier.
	CDFstatus status;            // CDF completion status.
	
	FILE*     fp;
		
	long cType;                 // Compression type
	long cParms[CDF_MAX_PARMS]; // Compression parameters
	
	status = CDFcreateCDF (fname.c_str(), &id);
	
	if (status != CDF_OK) 
		return false;
		
	long dims [ndims(M)];
	long dimv [ndims(M)];
	
	for (int i = 0; i < ndims(M); i++) {
		dims[i] = M.Dim(i);
		dimv[i] = VARY;
	}
	
	long imvn;
	
	status = CDFcreatezVar (id, dname.c_str(), CDF_FLOAT, 1L, (long) ndims(M), dims, VARY, dimv, &imvn);
	
	if (status != CDF_OK) 
		return false;
	
	cType = GZIP_COMPRESSION;
	cParms[0] = 5;             /* GZIP compression level */
	status = CDFsetzVarCompression (id, imvn, cType, cParms);
	
	if (status != CDF_OK) 
		return false;
	
	status = CDFcloseCDF (id);
	
	if (status != CDF_OK) 
		return false;
	
	return true;
	
#else 
	
	return false;
	
#endif
	
}


/**
 * @brief          Read matrix from CDF file.<br/> !!! FUNCTION NOT READY YET !!!
 * 
 * @param  M       Matrix
 * @param  fname   File name 
 * @param  dname   Data name
 * @param  dloc    Data location
 * @return         Success
 */
template <class T> inline static bool
CDFRead (Matrix<T>& M, const string fname, const string dname, const string dloc) {
	
	return false;
	
}


#endif //__IO_HPP__
