#ifdef PARC_MODULE_NAME

template <class T> long 
Matrix<T>::Import     (IceAs ias, long pos) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++)
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        _M[i+pos] = *data;
    
    return size;
    
}


template <class T> long 
Matrix<T>::Import(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs)")
        
    int i;
    
    for (i = 0; i < INVALID_DIM; i++)
        _dim[i] = (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    _M = new T[Size()]();
    nb_alloc = 1;
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < Size(); i++, data++)
        _M[i] = *data;
    
    return Size();
    
}


template <class T> long 
Matrix<T>::Export(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Export(IceAs)")
        
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); i++, data++)
        *data = _M[i];
    
    return Size();
    
}


template <class T> long
Matrix<T>::Export(IceAs ias, long pos) {

    ICE_SET_FN("Matrix<T>::Export(IceAs, long)")
        
        int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++) {
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    }
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        *data = _M[i+pos];
    
    return size;
    
}

#endif


