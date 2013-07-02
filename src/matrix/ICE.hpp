#ifndef MAX_ICE_DIM
#define MAX_ICE_DIM 16
#endif

size_t Import     (const IceAs& ias, const size_t pos) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < MAX_ICE_DIM; ++i)
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    T* data = (T*) ias.calcSplObjStartAddr();
    T* npos = &_M[pos];
    T* epos = npos+size; 
    
    for (; npos < epos; ++npos, ++data)
        *npos = *data;
    
    return size;
    
}


Matrix (const IceAs& ias) {
    ICE_SET_FN("Matrix<T,P>::Matrix(IceAs)")
    
        size_t n = MAX_ICE_DIM, i;
    _dim.resize(n);

    for (i = 0; i < MAX_ICE_DIM; ++i)
        _dim[i] = (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));

    // Remove trailing singleton dimensions
    for (i = 0; i < MAX_ICE_DIM; i++)
        if (_dim[i] == 1)
            n--;
        else
            n = MAX_ICE_DIM;
    
    // Resize skeleton
    _dim.resize(n);
    _res.resize(n,1.0);

    Allocate();

    T* data = (T*) ias.calcSplObjStartAddr() ;
    T* npos = &_M[0];
    T* epos = npos+Size(); 
    
    for (; npos < epos; ++npos, data++)
        *npos = *data;
    
}

size_t Import (const IceAs& ias) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs)")
        
    size_t n = MAX_ICE_DIM, i;

    _dim.resize(n);

    for (i = 0; i < MAX_ICE_DIM; ++i)
        _dim[i] = (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));

    // Remove trailing singleton dimensions
    for (i = 0; i < MAX_ICE_DIM; i++)
        if (_dim[i] == 1)
            n--;
        else
            n = MAX_ICE_DIM;
    
    // Resize skeleton
    _dim.resize(n);
    _res.resize(n,1.0);

    Allocate();

    T* data = (T*) ias.calcSplObjStartAddr() ;
    T* npos = &_M[0];
    T* epos = npos+Size(); 
    
    for (npos = 0; npos < epos; ++npos, data++)
        *npos = *data;
    
    return Size();
    
}


size_t Export (IceAs& ias) const {
    
    ICE_SET_FN("Matrix<T,P>::Export(IceAs)")
        
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); ++i, data++)
        *data = _M[i];
    
    return Size();
    
}


size_t Export (IceAs& ias, const size_t pos) const {

    ICE_SET_FN("Matrix<T,P>::Export(IceAs, long)")
        
        int  i    = 0;
    long size = 1;
    
    for (i = 0; i < MAX_ICE_DIM; ++i) {
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    }
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; ++i, data++)
        *data = _M[i+pos];
    
    return size;
    
}

