inline bool fexists (const char* fname) {
	
	std::ifstream fs (fname);
	
	if (fs)
		return true;
	else
		return false;

}

template<> inline std::ostream&  
Matrix<short>::Print     (std::ostream &os) const {

	for (int i = 0; i < _dim[LIN]; i++) {
		for(int j = 0; j < _dim[COL]; j++)
			printf ("%4i ", _M [i * _dim[COL] + j]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<double>::Print     (std::ostream &os) const {

	for (int i = 0; i < _dim[LIN]; i++) {
		for(int j = 0; j < _dim[COL]; j++)
			printf ("%+.4f ", _M [i * _dim[COL] + j]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<raw>::Print       (std::ostream& os) const {
	
	for (int i = 0; i < _dim[LIN]; i++) {
		for(int j = 0; j < _dim[COL]; j++)
			printf ("%+.4f+%+.4fi ", _M [i*_dim[COL]+j].real(), _M [i*_dim[COL]+j].imag());
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
bool Matrix<T>::dump (const char* fname) {
	
	int i = 0;

    if (fname != "") {
		
#ifdef HAVE_H5CPP_H

		try {
			
			Exception::dontPrint();
			
			hsize_t dims [INVALID_DIM];

			for (i = 0; i < INVALID_DIM; i++)
				dims[i] = _dim[INVALID_DIM-1-i];
			
			unsigned int mode = (fexists(fname)) ? H5F_ACC_RDWR : H5F_ACC_TRUNC;

			H5File    file  (fname, mode);

			if (typeid(T) == typeid(raw)) {

				DataSpace rspace (INVALID_DIM, dims);
				DataSpace ispace (INVALID_DIM, dims);
				FloatType type   (PredType::NATIVE_FLOAT);
				DataSet   iset = file.createDataSet(RAW_IM, type, ispace );
				DataSet   rset = file.createDataSet(RAW_RE, type, rspace );


				float     real[Size()];
				float     imag[Size()];
				
				for (i = 0; i < Size(); i++) {
					real[i] = raw(_M[i]).real();
					imag[i] = raw(_M[i]).imag();
				}
				

				rset.write  (real, type);
				iset.write  (imag, type);
				rset.close  ();
				iset.close  ();
				ispace.close();

			} else if (typeid(T) == typeid(double)) {

				DataSpace space (INVALID_DIM, dims);
				FloatType type  (PredType::NATIVE_DOUBLE);
				DataSet   set = file.createDataSet(HELPER, type, space);
				
				set.write  (_M, type);
				set.close  ();
				space.close();
			
			} else if (typeid(T) == typeid(short)) {

				DataSpace space (INVALID_DIM, dims);
				IntType   type  (PredType::NATIVE_SHORT);
				DataSet   set = file.createDataSet(PIXEL, type, space);
				
				set.write  (_M, type);
				set.close  ();
				space.close();
			
			}

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
		
		std::ofstream fout(fname , std::ios::out | std::ios::binary);
		
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
bool Matrix<T>::read (const char* fname) {
	
    if (fname != "") {
		
#ifdef HAVE_H5CPP_H

		try {
			
			Exception::dontPrint();
			
			H5File    file  (fname, H5F_ACC_RDONLY);
			DataSet dataset = file.openDataSet(RAW_RE);


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
		
        std::ifstream fin  (fname , std::ios::in | std::ios::binary);

        for (int i = 0; i < Size(); i++)
            fin.read  ((char*)(&(_M[i])), sizeof(T));
        
        fin.close();

#endif // HAVE_H5CPP_H

    }

	return true;
    
}

