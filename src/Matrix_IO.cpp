#include "Toolbox.hpp"
#include "tinyxml.h"
#include "mdhVB15.h"

#include <limits>
#include <map>


template <class T>
std::fstream& 
Matrix<T>::SkipLines (const std::ifstream& file, const unsigned int num) {
	
    file.seekg (std::ios::beg);
	
    for(int i=0; i < num - 1; ++i)
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	
    return file;
	
}


inline bool 
fexists (const char* fname) {
	
	std::ifstream fs (fname);
	
	if (fs)
		return true;
	else
		return false;
	
}


template<> inline std::ostream&  
Matrix<short>::Print     (std::ostream &os) const {
	
	for (int i = 0; i < _dim[COL]; i++) {
		for(int j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<double>::Print     (std::ostream &os) const {
	
	for (int i = 0; i < _dim[COL]; i++) {
		for(int j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<raw>::Print       (std::ostream& os) const {
	
	for (int i = 0; i < _dim[COL]; i++) {
		for(int j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f+%+.4fi ", _M [i + j * _dim[COL]].real(), _M [i + j * _dim[COL]].imag());
		printf("\n");
	}
	
	return os;
	
}


template <class T> std::ostream& 
operator<<    (std::ostream& os, Matrix<T>& M) {
	
	M.Print(os);
	return os;
	
}


template <class T>
bool Matrix<T>::pdump (std::string fname) {
	
	FILE *fp;
	
	if((fp=fopen(fname.c_str(), "wb"))==NULL) {
		printf("Cannot open %s file.\n", fname.c_str());
		return false;
	}
	
	if(fwrite(_M, sizeof(float), Size(), fp) != Size()) {
		printf("File read error.");
		return false;
	}
	
	fclose(fp);
	
	return true;

}

template <class T>
bool Matrix<T>::dump (std::string fname, std::string dname, std::string dloc, std::string fmt) {

	/*if (fmt.compare("mx"))
	  return mxdump (fname, dname, dloc);
	  else*/
	return h5dump (fname, dname, dloc);

}


template <class T>
bool Matrix<T>::h5dump (std::string fname, std::string dname, std::string dloc) {

	int i = 0;

	if (fname != "") {
		
#ifdef HAVE_H5CPP_H
		try {
			
#ifndef VERBOSE
			Exception::dontPrint();
#endif
			
			H5File        file; 
			
			try {
				file = H5File  (fname, H5F_ACC_TRUNC);
#ifdef VERBOSE
				printf ("File %s opened for RW\n", fname.c_str());
#endif
			} catch (Exception e) {
				file = H5File  (fname, H5F_ACC_TRUNC);
#ifdef VERBOSE
				printf ("File %s created for RW\n", fname.c_str());
#endif
			}
			
			Group group, *tmp;
			
			try {
				
				group = file.openGroup(dloc);
#ifdef VERBOSE
				printf ("Group %s opened for writing\n", dloc.c_str()) ;
#endif
				
			} catch (Exception e) {
				
				int    depth   = 0;

				std::vector<std::string> sv;

				Toolbox tb;
				tb.split(sv, dloc, "/");

				for (int i = 0; i < sv.size(); i++) {
					
					try {
						if (depth)
							group = (*tmp).openGroup(sv.at(i));
						else
							group = file.openGroup(sv.at(i));
					} catch (Exception e) {
						if (depth)
							group = (*tmp).createGroup(sv.at(i));
						else
							group = file.createGroup(sv.at(i));
					}
					
					tmp = &group;
					depth++;
					
				}
				
				group = (*tmp);
				
			}
			
			// One more field for complex numbers
			int tmpdim = (typeid(T) == typeid(raw)) ? INVALID_DIM+1 : INVALID_DIM;
			hsize_t* dims = new hsize_t[tmpdim];

			for (i = 0; i < INVALID_DIM; i++)
				dims[i] = _dim[INVALID_DIM-1-i];
			
			if (typeid(T) == typeid(raw))
				dims[INVALID_DIM] = 2;

			DataSpace space (tmpdim, dims);
			PredType*  type;
			
			delete [] dims;

			if (typeid(T) == typeid(raw)) {
				type = (PredType*) new FloatType (PredType::NATIVE_FLOAT);
				if (dname == "") 
					dname = "raw";
			} else if (typeid(T) == typeid(double)) {
				type = (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
				if (dname == "") 
					dname = "double";
			} else {
				type = (PredType*) new IntType   (PredType::NATIVE_SHORT);
				if (dname == "") 
					dname = "pixel";
			}
		
			DataSet set = group.createDataSet(dname, (*type), space);
				
			set.write   (_M, (*type));
			set.close   ();
			space.close ();
			file.close  ();
			
		} catch(FileIException      e) {
			e.printError();
			return false;
		} catch(DataSetIException   e) {
			e.printError();
			return false;
		} catch(DataSpaceIException e) {
			e.printError();
			return false;
		} catch(DataTypeIException  e) {
			e.printError();
			return false;
		}
		
#else // HAVE_H5CPP_H
		
		std::ofstream fout(fname.c_str() , std::ios::out | std::ios::binary);
		
		for (i = 0; i < Size(); i++) {
			std::cout << _M[i] << std::endl;
			fout.write ((char*)(&(_M[i])), sizeof(double));
		}
		
		fout.close();
		
#endif // HAVE_H5CPP_H

	}
	
	return true;
	
}



template <class T>
bool Matrix<T>::rsadjust (std::string fname) {

	int  dimk;
	long dimv;

	// XML file name, run cmd and repository
	std::stringstream    xmlf;
	std::stringstream    cmd;

	// XML objects
	TiXmlDocument*       doc = new TiXmlDocument();
	TiXmlElement*       meta;
	TiXmlNode*           fix;
	TiXmlNode*           vol;

	// Dimensions in XProtocol
	std::map < std::string, int > dims;
	dims["RawCol"] =  0; dims["RawLin"] =  1; dims["RawCha"] =  2; dims["RawSet"] =  3;
	dims["RawEco"] =  4; dims["RawPhs"] =  5; dims["RawRep"] =  6; dims["RawSeg"] =  7;
	dims["RawPar"] =  8; dims["RawSlc"] =  9; dims["RawIda"] = 10; dims["RawIdb"] = 11;
	dims["RawIdc"] = 12; dims["RawIdd"] = 13; dims["RawIde"] = 14; dims["RawAve"] = 15;

	// Create XML output from XProt
	cmd << "convert.pl ";
	cmd << fname;
	system (cmd.str().c_str());

	// Output filename
	xmlf << fname;
	xmlf << ".xml";

	// Read XML and go to rawobjectprovider for dimensions
	doc->LoadFile (xmlf.str().c_str());
	meta = doc->RootElement();
	vol  = meta->FirstChild ("XProtocol");
	vol  = vol->FirstChild ("ParamMap");
	vol  = vol->FirstChild ("ParamMap");
	fix  = vol->FirstChild ("Pipe");
	vol  = fix->FirstChild ("PipeService");
	vol  = fix->IterateChildren ("PipeService", vol);
	vol  = fix->IterateChildren ("PipeService", vol);
	fix  = fix->IterateChildren ("PipeService", vol);
	fix  = fix->FirstChild ("ParamFunctor");

	// Dimensions
	vol  = fix->FirstChild ("ParamLong");

	do {

		dimk = dims[((TiXmlElement*)vol)->Attribute( "name")];
		dimv = atoi(((TiXmlElement*)vol)->Attribute("value"));

		if (dimk == 0) dimv *= 2;
		if (dimv == 0) break;

		_dim[dimk] = dimv;

	} while ((vol = fix->IterateChildren ("ParamLong", vol))!=NULL);

	Reset();
	delete doc;

	return true;

}



template <class T>
bool Matrix<T>::rawread (std::string fname, std::string version) {
	
	FILE*         f;
	sMDH*         mdh;
	unsigned      l      = 0;
	unsigned long nscans = (Size() / _dim[COL]);
	unsigned      start;
	unsigned      size;

	// Get size information from XProtocol and resize 
	rsadjust(fname);

	// Assess data size
	f = fopen (fname.c_str(), "rb");
	fread (&l, sizeof(unsigned), 1, f);
	fseek (f,     l, SEEK_SET);
	start = ftell(f);
	fseek (f,    -1, SEEK_END);
	size  = ftell(f);
	fseek (f, start, SEEK_SET);
	
	//	printf ("Found %.1f MB of header and data.\n", (float) size / MB);
	printf ("Reading data from offset %i on ...\n", l );

	// Headers
	mdh  = (sMDH*) malloc (nscans * sizeof(sMDH));

	for (int i = 0; i < nscans; i++) {
		fread (&mdh[i], sizeof(sMDH), 1, f);
		fread (&_M[i * _dim[0]], sizeof (std::complex<float>), _dim[0], f);
		//if (i % 257 == 0)
		//printf ("%i\n", mdh[i].ushChannelId);
	}

	fclose (f);

	return true;


}

template <class T>
bool Matrix<T>::read (std::string fname, std::string dname, std::string dloc, std::string fmt) {

	//if (std::strcomp("mx", fmt.c_str()))
	//	return mxread (fname, dname, dloc);
	//else
		return h5read (fname, dname, dloc);
}


template <class T>
bool Matrix<T>::h5read (std::string fname, std::string dname, std::string dloc) {
	
    if (fname != "") {
		
#ifdef HAVE_H5CPP_H
		
#ifndef VERBOSE
		Exception::dontPrint();
#endif
		
		try {
			
			
			H5File    file (fname, H5F_ACC_RDONLY);
			DataSet   dataset = file.openDataSet(dloc+"/"+dname);
			DataSpace space   = dataset.getSpace();
			
			hsize_t*  dims    = (hsize_t*) malloc (space.getSimpleExtentNdims() * sizeof (hsize_t));
			int       ndim    = space.getSimpleExtentDims(dims, NULL);

			if (typeid(T) == typeid(raw)) 
				ndim--;

			for (int i = 0; i < ndim; i++)
				_dim[i] = dims[ndim-i-1];
			
			for (int i = ndim; i < INVALID_DIM; i++)
				_dim[i] = 1;
			
			if (nb_alloc)
				free (_M);

			_M = (T*) malloc (Size() * sizeof (T));
			

			std::cout << "rank: " << ndim << ", dimensions: ";
			for (int i = 0; i < ndim; i++) {
				std::cout << (unsigned long)(dims[i]);
				if (i == ndim - 1)
					std::cout << std::endl;
				else
					std::cout << " x ";
			}

			
			PredType*  type;
			
			if (typeid(T) == typeid(raw))
				type = (PredType*) new FloatType (PredType::NATIVE_FLOAT);
			else if (typeid(T) == typeid(double))
				type = (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
			else
				type = (PredType*) new IntType   (PredType::NATIVE_SHORT);
				
			dataset.read (_M, (*type));
			
			space.close();
			dataset.close();
			file.close();
			
		} catch(FileIException      e) {
			e.printError();
			return false;
		} catch(DataSetIException   e) {
			e.printError();
			return false;
		} catch(DataSpaceIException e) {
			e.printError();
			return false;
		} catch(DataTypeIException  e) {
			e.printError();
			return false;
		}
		
#else // HAVE_H5CPP_H
		
		std::ifstream fin  (fname.c_str() , std::ios::in | std::ios::binary);
		
		for (int i = 0; i < Size(); i++)
			fin.read  ((char*)(&(_M[i])), sizeof(T));
		
		fin.close();
		
#endif // HAVE_H5CPP_H
		
	}
	
	return true;
	
}

	
#ifdef HAVE_MAT_H

#include "mat.h"

template <class T>
bool Matrix<T>::mxread (std::string fname, std::string dname, std::string dloc) {

	
	// Open file ---------------------------------
	
	MATFile*  mf = matOpen (fname.c_str(), "r");
	
	if (mf == NULL) {
		printf ("Error opening file %s\n", fname.c_str());
		return false;
	}
	// -------------------------------------------
	
	// Get dimensions ----------------------------
	
	mxArray* mxa = matGetVariable (mf, dname.c_str());
	
	if (mxa == NULL) {
		printf ("Error opening variable %s\n", dname.c_str());
		return false;
	}
	
	mwSize        ndim = mxGetNumberOfDimensions(mxa);
	const mwSize*  dim = mxGetDimensions(mxa);
	
	for (int i = 0; i < ndim; i++)
		_dim[i] = (int)dim[i];
	
	for (int i = ndim; i < INVALID_DIM; i++)
		_dim[i] = 1;
	
	Reset();
	// -------------------------------------------
	
	// Copy from memory block ----------------------
	
	if (typeid(T) != typeid(short))
		memcpy(_M, mxGetPr(mxa), Size() * sizeof(T));
	else 
		for (int i = 0; i < Size(); i++)
			_M[i] = ((T*)mxGetPr(mxa))[i];
	// -------------------------------------------
	
	// Clean up and close file -------------------
	if (mxa != NULL)
		mxDestroyArray(mxa);
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
	// -------------------------------------------
	
}


template <class T>
bool Matrix<T>::mxdump (std::string fname, std::string dname, std::string dloc) {

	// Open file ---------------------------------

	MATFile*  mf = matOpen (fname.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", fname.c_str());
		return false;
	}
	// -------------------------------------------

	// Declare dimensions and allocate array -----

	mwSize   dim[INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++)
		dim[i] = (mwSize)_dim[i];
	
	mxArray*  mxa;

	if      (typeid(T) == typeid(double))
		mxa = mxCreateNumericArray (INVALID_DIM, dim, mxDOUBLE_CLASS,    mxREAL);
	else if (typeid(T) == typeid(raw))
		mxa = mxCreateNumericArray (INVALID_DIM, dim, mxDOUBLE_CLASS, mxCOMPLEX);
	else if (typeid(T) == typeid(short))
		mxa = mxCreateNumericArray (INVALID_DIM, dim,   mxINT8_CLASS,    mxREAL);
	// -------------------------------------------
	
	
	// Copy to memory block ----------------------
	
	if (typeid(T) != typeid(short))
		memcpy(mxGetPr(mxa), _M, Size() * sizeof(T));
	else
		for (int i = 0; i < Size(); i++)
			((T*)mxGetPr(mxa))[i] = _M[i];
	// -------------------------------------------
	
	// Write data --------------------------------
	int status = matPutVariable(mf, dname.c_str(), mxa);
	
	if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
		return false;
    }
	// -------------------------------------------
	
	
	// Clean up and close file -------------------
	
	if (mxa != NULL)
		mxDestroyArray(mxa);
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
	// -------------------------------------------

}

#endif
