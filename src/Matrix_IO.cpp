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
			printf ("%.4f ", _M [i * _dim[COL] + j]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<raw>::Print       (std::ostream& os) const {
	
	for (int i = 0; i < _dim[LIN]; i++) {
		for(int j = 0; j < _dim[COL]; j++)
			printf ("%.4f+%.4fi ", _M [i*_dim[COL]+j].real(), _M [i*_dim[COL]+j].imag());
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
	
    if (fname != "") {
		
#ifdef HAVE_H5CPP_H
		
		if (typeid(T) == typeid(raw))
			std::cout << raw(_M[0]).real() << std::endl;
		else 
			std::cout << _M[0] << std::endl;
		
		try {
			
			Exception::dontPrint();
			
			hsize_t dims [INVALID_DIM];

			for (int i = 0; i < INVALID_DIM; i++)
				dims[i] = _dim[i];
			

			H5File    file  (fname, H5F_ACC_TRUNC);

			if (typeid(T) == typeid(raw)) {

				DataSpace rspace (INVALID_DIM, dims);
				DataSpace ispace (INVALID_DIM, dims);
				FloatType type   (PredType::NATIVE_FLOAT);
				type.setOrder    (H5T_ORDER_LE);
				DataSet   iset = file.createDataSet( IMAG, type, ispace );
				DataSet   rset = file.createDataSet( REAL, type, rspace );


				float     real[Size()];
				float     imag[Size()];
				
				for (int j = 0; j < Size(); j++) {
					real[j] = raw(_M[j]).real();
					imag[j] = raw(_M[j]).imag();
				}
				

				rset.write (real, type);
				iset.write (imag, type);

				iset.close();
				ispace.close();

			} else {

				DataSpace rspace (INVALID_DIM, dims);
				
				FloatType type  (sizeof(T));
				type.setOrder   (H5T_ORDER_LE);
				
				DataSet   rset = file.createDataSet( REAL, type, rspace );
				
				rset.write (_M, type);
				rset.close();
				rspace.close();
			
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
		
		for (int i = 0; i < Size(); i++) {
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

        std::ifstream fin  (fname , std::ios::in | std::ios::binary);

        for (int i = 0; i < Size(); i++)
            fin.read  ((char*)(&(_M[i])), sizeof(T));
        
        fin.close();

    }

	return true;
    
}


