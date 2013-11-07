  /******************************
   ** definitions in oclMatrix **
   ******************************/


  /**
   * assignment operator
   */
  template <class T>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator=             (const oclMatrix <T> & mat)
  {
      
    // otherwise: self assignment
    if (this != &mat)
    {

      // copy members of Matrix <T>
      Matrix<T> :: _dim = mat._dim;
      Matrix<T> :: _res = mat._res;

      // temporary save current oclDataObject
      oclDataWrapper <T> * p_old_oclData = this -> mp_oclData;

      if (! mat.mp_oclData -> getLockState ())
      {

        /* if buffer is copied, CPU data need not to be copied */     /* TODO: if clause unnecessary */
        if (! mat.mp_oclData -> bufferCopyable ())
        {

          // copy data array of Matrix <T>
          Matrix<T> :: _M = mat.Container ();      // uses deep copy of valarray <T>/std::vector <T> (container)

        }
        else  /* but set size */
        {
          
          // copy data array of Matrix <T>
          Matrix<T> :: _M = mat.Container ();      // since resize also iterates over all elements !!!
          
        }
      
      
        // create new wrapper for gpu memory (and !new! cpu memory, copy object state!!!)
        /* important to copy object after CPU data (since memory pointer may become invalid) !!! */
        this -> mp_oclData = oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]), Matrix <T> :: Size (), * mat.mp_oclData, oclDataObject :: COPY_BUFFER);

      }
      else  /* -> so there are calculations on GPU  (bufferCopyable () should always return false) */
      {
      
        print_optional (class_name.c_str (), " :: operator=, KEEP_BUFFER", VERB_LOW);
      
        // create new wrapper for gpu memory (and !new! cpu memory, copy object state!!!)
        /* important to copy object after CPU data (since memory pointer may become invalid) !!! */
        this -> mp_oclData = oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]), Matrix <T> :: Size (), * mat.mp_oclData, oclDataObject :: KEEP_BUFFER);
      
      }
            
      // delete old oclDataObject
      delete p_old_oclData;
 
    }
        
    // for daisy chaining
    return *this;
       
  }




  /**
   * elementwise access (copy)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  operator[]              (const size_t & p)
  const
  {
  
    /* chose this order to avoid unnecessary data transfer from GPU */
    assert (p < Matrix <T> :: Size ());
      
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();                     /* TODO: possibility to copy just one element */
  
    }
      
    return Matrix <T> :: _M [p];
  
  }




  /**
   * elementwise access (reference)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  operator[]              (const size_t & p)
  {
  
    assert (p < Matrix <T> :: Size ());

    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {

      /* copy data to CPU */
      getData ();                     /* TODO: see above */

    }

    return Matrix <T> :: _M [p];
  
  }
  
  
  
  
  /**
   * data pointer (p-th position)
   */
  template <class T>
  inline
  const T *
  oclMatrix <T> ::
  Ptr                    (const size_t p)
  const
  {
  
    assert (p < this -> Size ());

    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {

      /* copy data to CPU */
      getData ();

    }
    
    return this -> _M.ptr (p);
  
  }




  /**
   * data object (lvalue)
   */
  template <class T>
  inline
  container <T> &
  oclMatrix <T> ::
  Container              ( )
  {

    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
    
    }
    
    return Matrix <T> :: _M;
  
  }
  
  
  
  
  /**
   * data object (rvalue)
   */
  template <class T>
  inline
  container <T>
  oclMatrix <T> ::
  Container                ( )
  const
  {
  
    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
    
    }
    
    return Matrix <T> :: _M;
  
  }
  
  
  
  
  /**
   * element at position p (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  At                      (const size_t & p)
  const
  {
  
    assert (p < this -> Size ());

    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
  
    }
    
    return Matrix <T> :: _M [p];
  
  }
  
  
  
  
  /**
   * element at position p (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  At                      (const size_t & p)
  {
  
    assert (p < this -> Size ());

    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
    
    }
    
    return Matrix <T> :: _M [p];
  
  }
  
  
  
  
  /**
   * element in slice (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  At                      (const size_t & x,
                           const size_t & y)
  const
  {
  
    assert (x < Matrix <T> :: _dim [0]);
    assert (y < Matrix <T> :: _dim [1]);
    
    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
  
    }
    
    return Matrix <T> :: _M [x + Matrix <T> :: _dim [COL] * y];
  
  }
  
  
  
  
  /**
   * element in slice (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  At                      (const size_t & x,
                           const size_t & y)
  {
  
    assert (x < Matrix <T> :: _dim [0]);
    assert (y < Matrix <T> :: _dim [1]);
    
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
  
    }

    return Matrix <T> :: _M [x + Matrix <T> :: _dim [COL] * y];
  
  }
  
  
  
  
  /**
   * element in volume (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  At                      (const size_t & x,
                           const size_t & y,
                           const size_t & z)
  const
  {
  
    assert (x < Matrix <T> :: _dim [0]);
    assert (y < Matrix <T> :: _dim [1]);
    assert (z < Matrix <T> :: _dim [2]);
    
    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
  
    }
    
    return Matrix <T> :: _M [x + Matrix <T> :: _dim [COL] * y];
  
  }
  
  
  
  
  /**
   * element in volume (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  At                      (const size_t & x,
                           const size_t & y,
                           const size_t & z)
  {
  
    assert (x < Matrix <T> :: _dim [0]);
    assert (y < Matrix <T> :: _dim [1]);
    assert (z < Matrix <T> :: _dim [2]);
    
    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
    
      /* copy data to CPU */
      getData ();
  
    }
    
    return Matrix <T> :: _M [x + Matrix <T> :: _dim [COL] * y];
  
  }
  
  
  
  
  /**
   * element in store (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  At                  (const size_t & col,
                       const size_t & lin,
                       const size_t & cha,
                       const size_t & set,
                       const size_t & eco,
                       const size_t & phs,
                       const size_t & rep,
                       const size_t & seg,
                       const size_t & par,
                       const size_t & slc,
                       const size_t & ida,
                       const size_t & idb,
                       const size_t & idc,
                       const size_t & idd,
                       const size_t & ide,
                       const size_t & ave)
  const
  {
  
    /* temporary pointer, for clear code */
    size_t * _dim = Matrix <T> :: _dim;
  
    assert (col < _dim [COL]);
    assert (lin < _dim [LIN]);
    assert (cha < _dim [CHA]);
    assert (set < _dim [SET]);
    assert (eco < _dim [ECO]);
    assert (phs < _dim [PHS]);
    assert (rep < _dim [REP]);
    assert (seg < _dim [SEG]);
    assert (par < _dim [PAR]);
    assert (slc < _dim [SLC]);
    assert (ida < _dim [IDA]);
    assert (idb < _dim [IDB]);
    assert (idc < _dim [IDC]);
    assert (idd < _dim [IDD]);
    assert (ide < _dim [IDE]);
    assert (ave < _dim [AVE]);
    
    /* copy data to CPU */
    getData ();
    
    return Matrix <T> :: _M [col+
                             lin*_dim[0]+
                             cha*_dim[0]*_dim[1]+
                             set*_dim[0]*_dim[1]*_dim[2]+
                             eco*_dim[0]*_dim[1]*_dim[2]*_dim[3]+
                             phs*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]+
                             rep*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]+
                             seg*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]+
                             par*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]+
                             slc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]+
                             ida*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]+
                             idb*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]+
                             idc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]+
                             idd*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]+
                             ide*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]+
                             ave*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]*_dim[14]
                            ];
  
  }    




  /**
   * element in store (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  At                  (const size_t & col,
                       const size_t & lin,
                       const size_t & cha,
                       const size_t & set,
                       const size_t & eco,
                       const size_t & phs,
                       const size_t & rep,
                       const size_t & seg,
                       const size_t & par,
                       const size_t & slc,
                       const size_t & ida,
                       const size_t & idb,
                       const size_t & idc,
                       const size_t & idd,
                       const size_t & ide,
                       const size_t & ave)
  {
  
    /* temporary pointer, for clear code */
    size_t * _dim = Matrix <T> :: _dim;

    assert (col < _dim [COL]);
    assert (lin < _dim [LIN]);
    assert (cha < _dim [CHA]);
    assert (set < _dim [SET]);
    assert (eco < _dim [ECO]);
    assert (phs < _dim [PHS]);
    assert (rep < _dim [REP]);
    assert (seg < _dim [SEG]);
    assert (par < _dim [PAR]);
    assert (slc < _dim [SLC]);
    assert (ida < _dim [IDA]);
    assert (idb < _dim [IDB]);
    assert (idc < _dim [IDC]);
    assert (idd < _dim [IDD]);
    assert (ide < _dim [IDE]);
    assert (ave < _dim [AVE]);
    
    /* copy data to CPU */
    getData ();
    
    return Matrix <T> :: _M [col+
                             lin*_dim[0]+
                             cha*_dim[0]*_dim[1]+
                             set*_dim[0]*_dim[1]*_dim[2]+
                             eco*_dim[0]*_dim[1]*_dim[2]*_dim[3]+
                             phs*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]+
                             rep*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]+
                             seg*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]+
                             par*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]+
                             slc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]+
                             ida*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]+
                             idb*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]+
                             idc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]+
                             idd*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]+
                             ide*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]+
                             ave*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]*_dim[14]
                            ];
  
  }




  /**
   * element in slice (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  operator()              (const size_t & x,
                           const size_t & y)
  const
  {
  
    return this -> At (x, y);
  
  }
  
  
  
  
  /**
   * element in slice (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  operator()              (const size_t & x,
                           const size_t & y)
  {
    
    return this -> At (x, y);
  
  }
  
  
  
  
  /**
   * element in volume (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  operator()              (const size_t & x,
                           const size_t & y,
                           const size_t & z)
  const
  {
    
    return this -> At (x, y, z);
  
  }
  
  
  
  
  /**
   * element in volume (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  operator()              (const size_t & x,
                           const size_t & y,
                           const size_t & z)
  {
  
    return this -> At (x, y, z);
  
  }
  
  
  
  
  /**
   * element in store (rvalue)
   */
  template <class T>
  inline
  T
  oclMatrix <T> ::
  operator()             (const size_t & col,
                          const size_t & lin,
                          const size_t & cha,
                          const size_t & set,
                          const size_t & eco,
                          const size_t & phs,
                          const size_t & rep,
                          const size_t & seg,
                          const size_t & par,
                          const size_t & slc,
                          const size_t & ida,
                          const size_t & idb,
                          const size_t & idc,
                          const size_t & idd,
                          const size_t & ide,
                          const size_t & ave)
  const
  {
  
    return this -> At (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
  
  }    




  /**
   * element in store (lvalue)
   */
  template <class T>
  inline
  T &
  oclMatrix <T> ::
  operator()            (const size_t & col,
                         const size_t & lin,
                         const size_t & cha,
                         const size_t & set,
                         const size_t & eco,
                         const size_t & phs,
                         const size_t & rep,
                         const size_t & seg,
                         const size_t & par,
                         const size_t & slc,
                         const size_t & ida,
                         const size_t & idb,
                         const size_t & idc,
                         const size_t & idd,
                         const size_t & ide,
                         const size_t & ave)
  {
  
    return this -> At (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
  
  }
  



  /**
   * unary plus
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator+               ()
  const
  {
    
    print_optional (class_name.c_str (), " :: operator+", op_v_level);
    
    return *this;
    
  }
  
  
  
  /**
   * unary minus
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator-               ()
  const
  {
  
    print_optional (class_name.c_str (), " :: operator-", op_v_level);
    
    // create matrix for result
    oclMatrix <T> inv = *this;
    
    // call operator function of oclOperations
    oclOperations <T, T> :: ocl_operator_mult_scalar (inv.mp_oclData, T (-1.), inv.mp_oclData, inv.Size ());
    
    return inv;
  
  }

  
  
 
  /**
   * elementwise addition of matrix
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator+               (const oclMatrix <S> & mat)
  const
  {

    print_optional (class_name.c_str (), " :: operator+", op_v_level);

    // create matrix for result
    oclMatrix <T> sum (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_add (this -> mp_oclData, mat.mp_oclData, sum.mp_oclData, this -> Size ());

    return sum;
    
  }



  /**
   * elementwise addition of scalar
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator+               (const S & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator+", op_v_level);
    
    // create matrix for result
    oclMatrix <T> sum = *this;
        
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_inc (sum.mp_oclData, s, sum.Size ());
    
    return sum;
  
  }




  /**
   * elementwise subtraction of matrix
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator-               (const oclMatrix <S> & mat)
  const
  {
    
    print_optional (class_name.c_str (), " :: operator-", op_v_level);
    
    // create matrix for result
    oclMatrix <T> diff (this -> Dim());
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_subtract (this -> mp_oclData, mat.mp_oclData, diff.mp_oclData, this -> Size ());
    
    return diff;
    
  }



  /**
   * elementwise subtraction of scalar
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator-               (const S & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator-", op_v_level);
    
    // create matrix for result
    oclMatrix <T> diff = *this;
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_dec (diff.mp_oclData, s, diff.Size ());
    
    return diff;
  
  }  
  
  
  
  /**
   * elementwise decrement by scalar
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator-=              (const S & dec)
  {
  
    print_optional (class_name.c_str (), " :: operator-=", op_v_level);
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_dec (this -> mp_oclData, dec, this -> Size ());
    
    return *this;
    
  }
  
  
  
  
  /**
   * elementwise decrement by matrix
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator-=              (const oclMatrix <S> & dec_mat)
  {
  
    print_optional (class_name.c_str (), " :: operator-=", op_v_level);
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_subtract (this -> mp_oclData, dec_mat.mp_oclData, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  

  /**
   * element increment by scalar
   */
  template <class T>  
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator+=              (const S & inc)
  {
  
    print_optional (class_name.c_str (), " :: operator+=", op_v_level);
 
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_inc (this -> mp_oclData, inc, this -> Size ());
  
    return *this;
  
  }
   
   
   
  
  /**
   * elementwise increment by matrix
   */
  template <class T>  
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator+=              (const oclMatrix <S> & inc_mat)
  {
  
    print_optional (class_name.c_str (), " :: operator+=", op_v_level);
 
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_add (this -> mp_oclData, inc_mat.mp_oclData, this -> mp_oclData, this -> Size ());
  
    return *this;
  
  }
  
  
  
  /**
   * elementwise raise to higher power
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator^               (const float & p)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator^", op_v_level);
    
    // create matrix for result
    oclMatrix <T> result (this -> Dim ());
        
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_pow (this -> mp_oclData, p, result.mp_oclData, this -> Size ());
        
    return result;
  
  }
  
  
  
  /**
   * elementwise raise to higher power and assignment
   */
  template <class T>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator^=              (const float & p)
  {
  
    print_optional (class_name.c_str (), " :: operator^=", op_v_level);
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_pow (this -> mp_oclData, p, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  
  /**
   * elementwise multiplication with a scalar
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator*               (const S & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator* (s)", op_v_level);
    
    // create result matrix
    oclMatrix <T> res (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_mult_scalar (this -> mp_oclData, s, res.mp_oclData, this -> Size ());
    
    return res;
  
  }
  
  
  
  /**
   * elementwise multiplication with a scalar
   * and assignment
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator*=              (const S & s)
  {
  
    print_optional (class_name.c_str (), " :: operator*= (s)", op_v_level);
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_mult_scalar (this -> mp_oclData, s, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  
  
  /**
   * elementwise multiplication with a matrix
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator*               (const oclMatrix <S> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator* (vector)", op_v_level);
    
    /* assert that matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create result matrix
    oclMatrix <T> res (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_mult_vector (this -> mp_oclData, mat.mp_oclData, res.mp_oclData, this -> Size ());
    
    return res;
  
  }
  
  
  
  /**
   * elementwise multiplication with a matrix
   * and assignment
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator*=              (const oclMatrix <S> & mat)
  {
  
    print_optional (class_name.c_str (), " :: operator*= (vector)", op_v_level);
    
    /* assert that matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // call operator function of oclOperations
    oclOperations <T, S> :: ocl_operator_mult_vector (this -> mp_oclData, mat.mp_oclData, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  
  /**
   * elementwise division by a scalar
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator/               (const S & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator/ (scalar)", op_v_level);
    
    /* create result matrix */
    oclMatrix <T> result (this -> Dim ());
    
    if (s != (S) 0)
      oclOperations <T, S> :: ocl_operator_div_scalar (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  /**
   * elementwise division by a scalar
   * and assignment
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator/=              (const S & s)
  {
  
    print_optional (class_name.c_str (), " :: operator/= (scalar)", op_v_level);
    
    if (s != (S) 0)
      oclOperations <T, S> :: ocl_operator_div_scalar (this -> mp_oclData, s, this -> mp_oclData, this -> Size ());
    else
      oclOperations <T> :: ocl_operator_assign (this -> mp_oclData, (T) 0, this -> Size ());
    
    return *this;
  
  }
  
  
  
  /**
   * elementwise division by a matrix
   */
  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator/               (const oclMatrix <S> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator/ (vector)", op_v_level);
    
    /* create result matrix */
    oclMatrix <T> result (this -> Dim ());
    
    oclOperations <T, S> :: ocl_operator_div_vector (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  
  /**
   * elementwise division by a matrix
   * and assignment
   */
  template <class T>
  template <class S>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator/=              (const oclMatrix <S> & mat)
  {
  
    print_optional (class_name.c_str (), " :: operator/= (vector)", op_v_level);
    
    oclOperations <T, S> :: ocl_operator_div_vector (this -> mp_oclData, mat.mp_oclData, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  
  
  /**
   * transposition of matrix
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator!               ()
  const
  {
  
    print_optional (class_name.c_str (), " :: operator*= (vector)", op_v_level);
    
    // load data to CPU
    getData ();
    
    // use CPU implementation of Matrix <T> for transposition
    return this -> Matrix <T> :: operator! ();
  
  }
  
  
  
  
  /**
   * scalar equality
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator==              (const T & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator==", op_v_level);

    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_equal (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
  
    return result;
  
  }



  /**
   * scalar inequality
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator!=          (const T & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator!=", op_v_level);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_inequal (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  /**
   * Scalar greater comparison.
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator>           (const T & s)
  const
  {

    print_optional (class_name.c_str (), " :: operator>", op_v_level);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_greater (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
    
  }


  /**
   * Scalar greater or equal comparison.
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator>=          (const T & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator>=", op_v_level);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_greater_equal (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
    
  }
      
      
  /**
   * Scalar less comparison.
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator<           (const T & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator<", op_v_level);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_less (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
    
  }
      
      
  /**
   * Scalar less or equal comparison.
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator<=           (const T & s)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator<=", op_v_level);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_less_equal (this -> mp_oclData, s, result.mp_oclData, this -> Size ());
    
    return result;
    
  }
  
  

  /**
   * elementwise equality with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator==          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator==", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_equal (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  /**
   * elementwise inequality with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator!=          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator!=", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_inequal (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }



  /**
   * elementwise greater with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator>           (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator>", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_greater (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }



  /**
   * elementwise greater or equal with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator>=          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator>=", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_greater_equal (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }



  /**
   * elementwise less with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator<           (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator<", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_less (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }



  /**
   * elementwise less or equal with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator<=          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator<=", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_less_equal (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  
  /**
   * masking (AND operation)
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator&           (const oclMatrix <cbool> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator&", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <T> result (this -> Dim ());
    
    oclOperations <T> :: ocl_operator_bitw_and (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  
  
  /**
   * elementwise AND operation with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator&&          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator&&", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_and (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  
  
  
  
  /**
   * elementwise OR operation with matrix
   */
  template <class T>
  oclMatrix <cbool>
  oclMatrix <T> ::
  operator||          (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator||", op_v_level);
    
    /* assert matrices have the same dimensions */
    assert (this -> _dim == mat._dim);
    
    // create matrix for result
    oclMatrix <cbool> result (this -> Dim ());
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_or (this -> mp_oclData, mat.mp_oclData, result.mp_oclData, this -> Size ());
    
    return result;
  
  }
  



  /**
   * matrix product
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  operator->*             (const oclMatrix <T> & fac_mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: operator->*", op_v_level);
    
    /* check for equality of inner dimensions */
    if (this -> _dim [1] != fac_mat._dim [0])
    {
      throw oclError (" Inner dimensions don't match!", "oclMatrix <T> :: operator->*");
    }
    
    oclMatrix <T> prod_mat (this -> _dim [0], fac_mat._dim [1]);
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_matprod (this -> mp_oclData, fac_mat.mp_oclData, prod_mat.mp_oclData,
                                                 this -> _dim [0],   this -> _dim [1],    fac_mat._dim [1],
                                                                0,                  0);
    
    return prod_mat;
  
  }

  

  /**
   * matrix product (transposed, right)
   */
  template <class T>  
  oclMatrix <T>
  oclMatrix <T> ::
  prodt               (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: prodt", op_v_level);
    
    return prod (mat, 'C');
  
  }



  /**
   * matrix product (transposed, right)
   */
  template <class T>  
  oclMatrix <T>
  oclMatrix <T> ::
  tprod               (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: tprod", op_v_level);
    
    return prod (mat, 'N', 'C');
  
  }
  
  
  
  /**
   * matrix product
   */
  template <class T>
  oclMatrix <T>
  oclMatrix <T> ::
  prod                (const oclMatrix <T> &    mat,
                       const      char     & transA,
                       const      char     & transB)
  const
  {
  
    print_optional (class_name.c_str (), " :: prod", op_v_level);
    
    assert (isxd (*this, 1) || is2d (*this));
    assert (isxd   (mat, 1) || is2d   (mat));
    
    int m, k, n,
        tr_A, tr_B;
    
    int cols_A = (int) this -> Dim (1), rows_A = (int) this -> Dim (0);
    int cols_B = (int)     mat.Dim (1), rows_B = (int)     mat.Dim (0);
    
    /* check for equality of inner dimensions */
    if      ( transA == 'N'                   &&  transB == 'N'                  ) assert (cols_A == rows_B);
    else if ( transA == 'N'                   && (transB == 'T' || transB == 'C')) assert (cols_A == cols_B);
    else if ((transA == 'T' || transA == 'C') &&  transB == 'N'                  ) assert (rows_A == rows_B);
    else if ((transA == 'T' || transA == 'C') && (transB == 'T' || transB == 'C')) assert (rows_A == cols_B);
    
    if (transA == 'N')
    {
      m = rows_A;
      k = cols_A;
    }
    else
    if (transA == 'T' || transA == 'C')
    {
      m = cols_A;
      k = rows_A;
    }

    if (transB == 'N')
      n = cols_B;
    else
    if (transB == 'T' || transB == 'C')
      n = rows_B;
    
    oclMatrix <T> prod_mat (m, n);
    
    // call operator function of oclOperations
    oclOperations <T> :: ocl_operator_matprod (this -> mp_oclData, mat.mp_oclData, prod_mat.mp_oclData,
                                                                m,              k,                   n,
                                                                0,              0);
    
    return prod_mat;
  
  }
  
  
  
  /**
   * Dot product.
   */
  template <class T>
  T
  oclMatrix <T> ::
  dot                 (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: dot", op_v_level);
    
    /* check dimensions */
    assert (this -> Dim (1) == 1);
    assert (mat.Dim (1) == 1);
    assert (this-> Dim (0) == mat.Dim (0));
    
    /* call matrix matrix product */
    oclMatrix <T> prod_mat = prod (mat, 'T');
    
    prod_mat.getData ();
    return * (prod_mat.Ptr ());
  
  }



  /**
   * Dot product.
   */
  template <class T>
  T
  oclMatrix <T> ::
  dotc                (const oclMatrix <T> & mat)
  const
  {
  
    print_optional (class_name.c_str (), " :: dotc", op_v_level);
    
    /* check dimensions */
    assert (this -> Dim (1) == 1);
    assert (mat.Dim (1) == 1);
    assert (this-> Dim (0) == mat.Dim (0));
    
    /* call matrix matrix product */
    oclMatrix <T> prod_mat = prod (mat, 'C');
    
    prod_mat.getData ();
    return * (prod_mat.Ptr ());
  
  }



  /**
   * cast operator
   */
  template <class T>
  template <class S>
  oclMatrix <T> ::
  operator oclMatrix <S> ()
  const
  {
  
    print_optional (class_name.c_str (), " :: operator oclMatrix <S>", op_v_level);
    
    oclMatrix <S> cast (this -> Dim ());
    
    oclOperations <T, S> :: ocl_operator_cast (this -> mp_oclData, cast.mp_oclData, this -> Size ());
    
    return cast;
  
  }
  




  /**
   * resize matrix to  m x n.
   */
  template <class T>
  inline
  void
  oclMatrix <T> ::
  Resize              (const size_t & m,
                       const size_t & n)
  {
  
    print_optional (class_name.c_str (), " :: Resize", op_v_level);
  
    /* need to update CPU data ?? */
    if (this -> mp_oclData -> getMemState ()
     && ! this -> mp_oclData -> getCPUModified ())
    {
  
      /* copy data to CPU */
      getData ();

    }
    
    /* function of base class */
    Matrix <T> :: Resize (m, n);

    /* delete old data object */
    delete this -> mp_oclData;
    
    /* create new data object */
    this -> mp_oclData = oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]), Matrix <T> :: Size ());
  
  }
  
  
  
  
  /**
   * Purge data and free RAM.
   */
  template <class T>
  inline
  void
  oclMatrix <T> ::
  Clear               ( )
  {
  
    print_optional (class_name.c_str (), " :: Clear", op_v_level);
  
    /* call function of base class */
    Matrix <T> :: Clear ( );
    
    /* delete data object */
    delete this -> mp_oclData;
  
  }
  
  
  
  /**
   * assignment operator for valarray
   */
  template <class T>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator=           (const std::valarray <T> & v)
  {
  
    print_optional (class_name.c_str (), " :: operator= (valarray)", op_v_level);
  
    /* valarrays must have the same size */
    assert (Matrix <T> :: _M.size () == v.size ());
    
    /* check if valarrays are the same */
    if (& Matrix <T> :: _M != & v)
    {
    
      /* check if data on GPU are up to date */
      if (this -> mp_oclData -> getMemState ()
        && ! this -> mp_oclData -> getCPUModified ())
      {
    
        // wait for calculations to end
        while (this -> mp_oclData -> getLockState ())
          usleep (5);
      
        // create data object
        oclDataObject * cp_obj = oclOperations <T> :: make_GPU_Obj (&(v [0]), v.size ());
    
        // launch kernel operator
        oclOperations <T> :: ocl_operator_copy (this -> mp_oclData, cp_obj, this -> Size ());

      }
      else
      {
    
        /* copy on CPU */
        Matrix <T> :: _M = v;
    
      }
  
    }
    
    return *this;
  
  }
      
      
      
      
  /**
   * assignment operator for scalar
   */
  template <class T>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator=           (const T & s)
  {
  
    print_optional (class_name.c_str (), " :: operator= (scalar)", op_v_level);
  
    /* check if data on GPU are up to date */
    if (this -> mp_oclData -> getMemState ()
      && ! this -> mp_oclData -> getCPUModified ())
    {
    
      // wait for calculations to end
      while (this -> mp_oclData -> getLockState ())
        usleep (5);
    
      // launch kernel operator
      oclOperations <T> :: ocl_operator_assign (this -> mp_oclData, s, this -> Size ());

    }
    else
    {
    
      /* copy on CPU */
      Matrix <T> :: _M = s;
      
    }
    
    return *this;
    
  }
      
      
      
  
  /**
   * copy data to CPU, if
   *  - GPU buffer exists
   *  - data are not synchronized
   */
  template <class T>
  inline
  void
  oclMatrix <T> ::
  getData                 ()
  const
  {
  
    print_optional (class_name.c_str (), " :: getData", VERB_LOW);
  
    // check if data on CPU are up to date
    if (! this -> mp_oclData -> getSyncState ()
       && this -> mp_oclData -> getMemState ())
    {
  
      print_optional (" getData for object: ", mp_oclData -> getID (), VERB_LOW);
  
      int i = 0;
      bool done = false;
      while (! done)
      {
      
        try
        {
        
          /* print objects state */
//          this -> mp_oclData -> print ();
        
          // load GPU data to CPU
          this -> mp_oclData -> getData ();
          
          done = true;

        }
        catch (oclError & err)
        {
        
          std::cout << " wait: " << i++ << "\n " << err << std::endl;
        
        }
      
      }

    }
  
  }



template <class T>
oclError &
assign                    (          Matrix <T> &     mat,
                            const oclMatrix <T> & ocl_mat )
{

  print_optional (" :: assign", VERB_HIGH);

  try
  {

    /* copy data to CPU */
    ocl_mat.getData ();
  
    /* call assignment operator of base class */
    mat = ocl_mat;

  }
  catch (oclError & err)
  {
  
    throw oclError (err, ":: assign (Matrix <T>, oclMatrix <T>)");
  
  }

}



