#ifndef __MATRIX_H__
#define __MATRIX_H__

#ifdef ICEIDEAFUNCTORS_EXPORTS

  #include     "MrServers/MrVista/include/Ice/IceBasic/IceAs.h"
  #include     "MrServers/MrVista/include/Ice/IceBasic/IceObj.h"
  #include     "MrServers/MrVista/include/Parc/Trace/IceTrace.h"

#else

enum IceDim {
    COL, LIN, CHA, SET, ECO, PHS, REP, SEG, PAR, SLC, IDA, IDB, IDC, IDD, IDE, AVE, INVALID_DIM
};

#endif

#include <complex>
#include <iostream>
#include <fstream>
#include <assert.h>


/**
 * Short test if the matrix is a vector.
 */
# define VECT(M) assert((M)->width() == 1 || (M)->height() == 1);

/**
 * Return absolute value.
 */
# define ABS(A) (A > 0 ? A : -A)

/**
 * Return rounded value as in MATLAB.
 * i.e. ROUND( 0.1) ->  0
 *      ROUND(-0.5) -> -1
 *      ROUND( 0.5) ->  1
 */
# define ROUND(A) ( floor(A) + ((A - floor(A) >= 0.5) ? (A>0 ? 1 : 0) : 0))

/**
 * Defined for memory allocation testing.
 */
#define ALLOC
#ifdef ALLOC
  static int nb_alloc = 0;
#endif

/**
 * Old friend pi.
 */
# define PI  3.1415926535897931159979634685441851615906

/**
 * @brief   Matrix template
 * 
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */

template <typename T> class Matrix {
    

public:
    
    
    /**
     * @name Constructors and destructors
     *       Constructors and destructors
     */
    //@{
    

    /**
     * @brief           Default constructor.
     */
    inline              
    Matrix              ();
    
    
    /**
     * @brief           Constructs a new 1^16 matrix initialized with scalar s.
     *
     * @param  s        Scalar Value.
     * @return          A new matrix with one cell.
     */
    inline              
    Matrix              (T s) {};
    
    
    /**
     * @brief           Constructs an uninitialised matrix with desired size.
     *
     * @param  col      Columns
     * @param  lin      Lines
     * @param  cha      Channels
     * @param  set      Sets
     * @param  eco      Echoes
     * @param  phs      Phases
     * @param  rep      Repetitions
     * @param  seg      Segments
     * @param  par      Partitions
     * @param  slc      Slices
     * @param  ida      IDA
     * @param  idb      IDB
     * @param  idc      IDC
     * @param  idd      IDD
     * @param  ide      IDE
     * @param  ave      Averages
     * @return          Constructed matrix.
     */
    inline              
    Matrix               (int col, int lin, int cha, int set, 
                          int eco, int phs, int rep, int seg, 
                          int par, int slc, int ida, int idb, 
                          int idc, int idd, int ide, int ave);
    
    
    /**
     * @brief           Constructs an uninitialised matrix with desired size in integer matrix.
     *
     * @param  dim      integer array of 16 elements with dimensions in following order:
     *                  COL, LIN, CHA, SET, ECO, PHS, REP, SEG, PAR, SLC, IDA, IDB, IDC, IDD, IDE, AVE
     * @return          Constructed matrix.
     */
    inline              
    Matrix              (Matrix<int> &dim);
    

    #ifdef ICEIDEAFUNCTORS_EXPORTS
    
    /**
     * @brief           Reset and fill data from IceAs
     *                   
     * @param  ias      IceAs containing data
     * 
     * @return          Amount of data read
     */
    inline long         
    Import              (IceAs ias);

    /**
     * @brief           Import an amount of data from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @param  pos      Import data starting at position pos own repository
     * 
     * @return          Amount of data read
     */
    inline long         
    Import              (IceAs ias, long pos);

    /**
     * @brief           Export data to ias
     *                   
     * @param  ias      IceAs for data export
     * 
     * @return          Amount of data exported
     */
    inline long         
    Export              (IceAs ias);

    /**
     * @brief           Partially export data to ias 
     * 
     * @param  ias      IceAs for data export
     * @param  pos      Export data starting at position pos of our repository
     */
    inline long         
    Export              (IceAs ias, long pos);
 
    #endif


    /**
     * @brief           Default destructor.
     */
    ~Matrix             ();
    
    //@}


    /**
     * @name            Elementwise access
     *                  Elementwise access
     */
    
    //@{
    
    
    /**
     * @brief           Get an element of a vector at position p, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested Scalar value.
     */
    T                   
    operator[]          (int p)                             const ;
    
    
    /**
     * @brief           Get an element of a vector at position p, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Reference to the requested Scalar value.
     */
    T                   
    &operator[]         (int p)                              ;

    
    /**
     * @brief           Get value in slice
     *  
     * @param  col      Column
     * @param  lin      Line
     *
     * @return          Value
     */
    inline T            
    at                  (int col, int lin)  const {
        return _M[col + _dim[LIN]*lin ];
    }

    
    /**
     * @brief            Set value in slice
     *  
     * @param  col       Column
     * @param  lin       Line
     *
     * @return           Reference to element
     */
    inline T&           
    at                  (int col, int lin) {
        return _M[col + _dim[LIN]*lin ];
    }

    
    /**
     * @brief            Get value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Value
     */
    inline T            
    at                   (int col, int lin, int slc)  const {
        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    
    /**
     * @brief            Set value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Reference to element
     */
    inline T&            
    at                   (int col, int lin, int slc) {
        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    
    /**
     * @brief           Get the element at position p of the vector, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T                  
    operator()          (int p) const;

    
    /**
     * @brief           Get the element at position p of the vector, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (int p) ;

    
    /**
     * @brief            Get value in slice
     *
     * @param  col       Column
     * @param  lin       Line
     *
     * @return           Value
     */
    T                  
    operator()           (int col, int lin) const {
        return _M[col + _dim[LIN]*lin ];
    }
    

    /**
     * @brief            Set value in slice
     *
     * @param  col       Column
     * @param  lin       Line
     *
     * @return           Reference to element.
     */
    T&                  
    operator()           (int col, int lin) {
        return _M[col + _dim[LIN]*lin ];
    }
    
    
    /**
     * @brief            Get value in volume
     *
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Value
     */
    T                  
    operator()           (int col, int lin, int slc) const {
        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    
    /**
     * @brief            Set value in column
     *
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Reference to element
     */
    T&                 
    operator()           (int col, int lin, int slc) {
           return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    };
    
    //@}
    




    /**
     * @name            Partial access functions.
     *                  Functions for access to parts of data.
     */
    
    //@{
    
    /**
     * @brief           Get a vector with as copy of line l of this Matrix, i.e. this(l,:).
     *
     * @param  l        Requested line.
     * @return          Vector copied from requested row.
     */
    Matrix<T>          
    lin                 (int l)                             const;
    
    
    /**
     * @brief           Get a vector with as copy of column c of this Matrix, i.e. this(:,c).
     *
     * @param  c        Requested column.
     * @return          Vector copied from requested column.
     */
    Matrix<T>           
    col                 (int c)                             const;
    
    
    /**
     * @brief           Get a matrix as a copy of lines l[i] of this matrix, i.e. this(l,:);
     *
     * @param  l        Requested rows.
     * @return          Matrix copied from requested rows.
     */
    Matrix<T>           
    lin                 (Matrix<int> l)                     const;
    
    
    /**
     * @brief           Get a matrix as a copy of columns c[i] of this matrix, i.e. this(:,c);
     *
     * @param  c        Requested rows.
     * @return          Matrix copied from requested columns.
     */
    Matrix<T>           
    col                 (Matrix<int> c)                     const;
    
    
    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    int                 
    height              () const {return _dim[LIN];}
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    int                 
    width               () const {return _dim[COL];}
    
    
    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    int                 
    m                   () const {return _dim[LIN];}


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    int                 
    n                   () const {return _dim[COL];}


    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline int          
    Dim                 (int i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline const int*   
    Dim                 ()                                     const {return _dim;}
    
    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline void         
    Dim                 (const int* dim)                                     const {
        for (int i = 0; i<INVALID_DIM; i++)
            _dim[i] = dim[i];
    }
    
    
    /**
     * @brief           Resize to dim and zero
     *
     * @param  dim      New dimensions
     */
    inline void         
    Reset               (const int* dim)                                      {
        for (int i = 0; i < INVALID_DIM; i++)
            _dim[i] = dim[i];
        //if (_M != 0)
        //delete [] (_M);
        _M = new T[Size()]();
    }
    

    /**
     * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
     *
     * @return          Number of cells.
     */
    long                
    Size                ()                                    const;
    
    
    /**
     * @brief           Get the number of matrix cells, i.e. Size * sizeof(T).
     *
     * @return          Size in RAM in bytes.
     */
    int                 
    SizeInRAM           ()                                    const;

    //@}
    
    

    /**
     * @name            Operators
     *                  Operator definitions. Needs big expansion still.
     */
    
    //@{
    

    /**
     * @brief           Get one or more elements of a vector, i.e. this(p).
     *
     * @param  p        Requested positions
     * @return          Vector of requested elements.
     */
    Matrix<T>           
    operator()          (Matrix<int> &p)                    const;

    
    /**
     * @brief           Get one or more elements of a matri, i.e. this=[1 2 3;4 5 6;7 8 9],
     *                  m1=[2], m2=[1 3], returns [2 8].
     *
     * @param  r        Requested rows.
     * @param  c        Requested columns.
     * @return          matrix of requested elements.
     */
    Matrix<T>           
    operator()          (Matrix<int> &c, Matrix<int> &l)  const;
    
    
    /**
     * @brief           Assignment operator. i.e. this = m.
     *
     * @param  M        The assigned matrix.
     */
    Matrix<T>           
    operator=           (Matrix<T> &M);
    
    
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator->*         (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise substruction of two matrices. i.e. this - M.
     *
     * @param  M        The reference to the substruent.
     */
    Matrix<T>           
    operator-           (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise substruction all elements by a given value. i.e. this - s.
     */
    Matrix<T>           
    operator-           (T s);
    
    
    /**
     * @brief           Substruct from 0. i.e. 0 - this.
     */
    Matrix<T>           
    operator-           ();
    
    
    /**
     * @brief           Transposition. i.e. this'.
     */
    Matrix<T>           
    operator!           () const;
    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     */
    Matrix<T>           
    operator&           (Matrix<bool> &M);
    
    
    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          True if and only if all elements of the matrix equal s.
     */
    Matrix<bool>        
    operator==          (T s);
    
    
    /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator!=          (T s);
    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator>           (T s);
    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator>=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator<=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator<           (T s);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator==          (Matrix<T> M);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator!=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator>=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     */
    Matrix<bool>        
    operator<=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator>           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator<           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, tel que result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<T>           
    operator||          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, tel que result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<T>           
    operator&&          (Matrix<T> &M);


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Exponent.
     */
    Matrix<T>           
    operator^           (int p);
    
    //@}





    /**
     * @name Data manipulation functions.
     *       Data manipulation functions.
     */
    
    //@{
    
    /**
     * @brief           Matrix multiplication. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator*           (Matrix<T> &M);
    
    /**
     *                  Elementwise multiplication with a scalar. i.e. this * m.
     */
    Matrix<T>           
    operator*           (T s);
    
    //@}
    
    
    
    
    /**
     * @name Display and IO functions.
     *       Display and IO functions.
     */
    
    //@{
    

    /**
     * @brief           Print contents to output stream.
     *
     * @param  os       The output stream.
     * @return          The output stream.
     */
    std::ostream        
    &print              (std::ostream &os) const;
    

    /**
     * @brief           Print contents to output stream.
     *
     * @param  os       The output stream.
     * @return          The output stream.
     */
    std::ostream        
    &operator<<         (std::ostream &os);


    /*
     * @brief           Dump binary column-major matrix.
     * 
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    dump                (char* fname);
    

    /**
     * @brief           Read in binary column-major matrix.
     *
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    read                (char* fname);
    
    //@}
    
    
    
    /**
     * @name            Other functions.
     *                  Other functions.
     */
    
    //@{
    
    
    /**
     * @brief           Greatest element of the matrix. i.e. max(M0).
     *
     * @return          The scalar value of the greatest element.
     */
    T                   
    Max                 ();    
    
    
    /**
     * @brief           Smallest element of the matrix. i.e. min(M0).
     *
     * @return          The scalar value of the smallest element.
     */
    T                   
    Min                 ();
    
    
    /**
     * @brief           Get maximum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Maxabs              ();
    
    
    /**
     * @brief           Get minimum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Minabs              ();
    
    /**
     * @brief           Matrix Product.
     *
     * @return          Product of this and M.
     */
    Matrix<T>           
    prod                (Matrix<T> &M);
    
    /**
     * @brief           Transposition.
     *
     * @return          The transposed matrix
     */
    Matrix<T>           
    tr()                const;
    
    //@}
    
    

private:
    
    int                 _dim[INVALID_DIM]; /// Dimnesions
    T*                  _M;                /// Data repository
    
};


template <class T> 
Matrix<T>::Matrix () {

    for (int i = 0; i < INVALID_DIM; i++)
        _dim [i] = 1;

    #ifdef ALLOC
        nb_alloc++;
    #endif

};


template <class T> 
Matrix<T>::Matrix (int col, int lin, int cha, int set, 
                   int eco, int phs, int rep, int seg, 
                   int par, int slc, int ida, int idb, 
                   int idc, int idd, int ide, int ave) {

    _dim[COL] = col;
    _dim[LIN] = lin;
    _dim[CHA] = cha;
    _dim[SET] = set;
    _dim[ECO] = eco;
    _dim[PHS] = phs;
    _dim[REP] = rep;
    _dim[SEG] = seg;
    _dim[PAR] = par;
    _dim[SLC] = slc;
    _dim[IDA] = ida;
    _dim[IDB] = idb;
    _dim[IDC] = idc;
    _dim[IDD] = idd;
    _dim[IDE] = ide;
    _dim[AVE] = ave;

    #ifdef ALLOC
        nb_alloc++;
    #endif

    _M = new T[Size()]();

};


template <class T> 
Matrix<T>::Matrix (Matrix<int> &dim) {

    assert (dim.Size() == INVALID_DIM);

    for (int i = 0; i < INVALID_DIM; i++) {
        _dim[i] = dim[i];
    }
                             
    #ifdef ALLOC
        nb_alloc++;
    #endif

    _M = new T[Size()]();

};


#ifdef ICEIDEAFUNCTORS_EXPORTS

template <class T>
long Matrix<T>::Import(IceAs ias, long pos) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++)
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        _M[i+pos] = *data;
    
    return size;
    
};

template <class T> 
long Matrix<T>::Import(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs)")
        
        int i;
    
    for (i = 0; i < INVALID_DIM; i++)
        _dim[i] = (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    #ifdef ALLOC
        nb_alloc++;
    #endif
    _M = new T[Size()]();
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < Size(); i++, data++)
        _M[i] = *data;
    
    return Size();
    
};

template <class T> 
long Matrix<T>::Export(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Export(IceAs)")
        
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); i++, data++)
        *data = _M[i];
    
    return Size();
    
};

template <class T> 
long Matrix<T>::Export(IceAs ias, long pos) {

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
    
};



#endif

template <class T> 
Matrix<T>::~Matrix() {
    
#ifdef ICEIDEAFUNCTORS_EXPORTS
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

    if (_M !=0 )
    	delete [] (_M);
    
    #ifdef ALLOC
        nb_alloc--;
    #endif
    
};


template <class T>
Matrix<T> Matrix<T>::operator* (Matrix<T> &M) {
    
    assert(Size() == M.Size());
    
    Matrix<T> res(_dim);
    
    for (int i = 0; i < Size(); i++)
        res[i] = _M[i] * M[i];
    
    return res;
    
};


template <class T>
long         Matrix<T>::Size() const {
    
    long size = 1;
    
    for (int i = 0; i < INVALID_DIM; i++)
        size *= _dim[i];
    
    return size;
    
};


template <class T>
inline int  Matrix<T>::SizeInRAM() const {
    
    return Size() * sizeof(T);
    
};


template <class T>
T           Matrix<T>::operator[]  (int p) const {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
};

template <class T>
T           &Matrix<T>::operator[] (int p) {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
};


template <class T>
Matrix<T>   Matrix<T>::operator()(Matrix<int> &m) const {
    
    VECT(&m);
    
    Matrix<T> res(1, m.Size());
    
    for (int i = 0; i < res.Size(); i++)
        res[i] = (*this)[m[i]];
    
    return res;

};


template <class T>
Matrix<T> Matrix<T>::operator()(Matrix<int> &m1, Matrix<int> &m2) const {

    Matrix<T> res(m1.Size(),  m2.Size());
    
    for (int i = 0; i < m1.Size(); i++)
        for (int j = 0; j < m2.Size(); j++)
            res[i * m2.Size() + j] = _M[m1[i] * _dim[LIN] + m2[j]];
    
    return res;

};


template <class T>
T Matrix<T>::operator() (int a) const {

    assert(a >= 0);
    assert(a <  Size());

    return _M[a];

};


template <class T>
T Matrix<T>::Max() {

    T old = _M[0];

    for (int i = 0; i < Size(); i++)
        if (_M[i] > old)
            old = _M[i];

    return old;

}


template <class T>
T  Matrix<T>::Maxabs() {

    T old = fabs(_M[0]);

    for (int i = 0; i < Size(); i++)
        if (fabs(_M[i]) > old)
            old = fabs(_M[i]);

    return old;

}


template <class T>
T Matrix<T>::Min() {

    T old = _M[0];

    for (int i = 0; i < Size(); i++)
        if (_M[i] < old)
            old = _M[i];

    return old;

}


template <class T>
T  Matrix<T>::Minabs() {

    T old = fabs(_M[0]);

    for (int i = 0; i < Size(); i++)
        if (fabs(_M[i]) < old)
            old = fabs(_M[i]);

    return old;

}



template <class T>
Matrix<T> Matrix<T>::operator=(Matrix<T> &M) {
    
    int i;

    for (i = 0; i < INVALID_DIM; i++) 
        _dim[i] = M.Dim()[i];
    
    //if (_M != 0)
        //delete[](_M);

    _M = new T[Size()]();

    for (i = 0; i < Size(); i++)
        _M[i] = M[i];

    return *this;

}


template <class T>
Matrix<bool> Matrix<T>::operator==(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>=(T s) {

    Matrix<bool> res(_dim);
    
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    
    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator!=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] > s);

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator->*(Matrix<T> &M) {

    return this->prod(M);

}

template <class T>
Matrix<T> Matrix<T>::operator!() const {

    return this->tr();

}

template <class T>
Matrix<T> Matrix<T>::operator&(Matrix<bool> &M)    {

    for (int i = 0; i < INVALID_DIM; i++) 
        assert (_dim[i] == (int) M.Dim()[i]);


    int k = 0, i;
    for (i = 0; i < Size(); i++)
        if (M[i])
            k++;
    
    Matrix<T> res(k, 1);
    
    k = 0;
    for (i = 0; i < Size(); i++)
        if (M[i]) {
            res[k] = i;
            k++;
        }

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator&&(Matrix<T> &M) {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator||(Matrix<T> M) {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (M[i] || _M[i]);

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator==(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator>=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator<=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator!=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator>(Matrix<T> M) {
    assert(_dim[COL] == M.wridth() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator<(Matrix<T> M) {

    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-() {

    Matrix<T> res(_dim[LIN], _dim[COL]);

    for (int i = 0; i < Size(); i++)
        res[i] = -_M[i];

    return res;

};

template <class T>
static T power(T p, int b) {

    assert(b > 0);

    T res = p;

    b--;

    while (b > 0) {
        res *= p;
        b--;
    }

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator^(int p) {
    Matrix<T> res(_dim[LIN], _dim[COL]);
    for (int i = 0; i < Size(); i++)
        if (p == 0)
            res[i] = 1;
        else
            if (p > 0)
                res[i] = power(_M[i], p);
            else
                res[i] = 1 / power(_M[i], -p);
    return res;
}


template <class T>
Matrix<T> Matrix<T>::prod(Matrix<T> &m) {
    
    assert(_dim[COL] == m.height());
    Matrix<T> res(_dim[LIN], m.width());

    for (int i = 0; i < res.height(); i++)
        for (int j = 0; j < res.width(); j++) {
            res[i * res.width() + j] = 0;
            for (int k = 0; k < _dim[COL]; k++)
                    res[i * res.width() + j] += _M[i * _dim[COL] + k] *
                    m[k * m.width() + j];
        }
    return res;
}

template <class T>
Matrix<T> Matrix<T>::tr() const {

    Matrix<T> res(_dim);

    for (int i = 0; i < res.height(); i++)
        for (int j = 0; j < res.width(); j++)
            res[i * res.width() + j] = _M[j * _dim[COL] + i];

    return res;

}

template <class T>
bool Matrix<T>::dump (char* fname) {

	// TODO: Error checking.
	//       File not found.
	//       HDF5 file format.

    if (fname != "") {

        std::ofstream fout(fname , std::ios::out | std::ios::binary);

        for (int i = 0; i < Size(); i++)
            fout.write ((char*)(&(_M[i])), sizeof(T));
        
        fout.close();

    }
    
}

template <class T>
bool Matrix<T>::read (char* fname) {

	// TODO: Error checking.
	//       File not found.
	//       EOF before Size.
	//       HDF5 file format.

    if (fname != "") {

        std::ifstream fin  (fname , std::ios::in | std::ios::binary);

        for (int i = 0; i < Size(); i++)
            fin.read  ((char*)(&(_M[i])), sizeof(T));
        
        fin.close();

    }
    
}

#endif // __MATRIX_H__


