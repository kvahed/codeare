# ifndef __OCL_MATRIX_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_MATRIX_HPP__




  /**************
   ** includes **
   **************/
 
  // CoDEARE
  # include "../../Matrix.hpp"
  
  // ocl
  # include "../oclSettings.hpp"
  # include "../oclConnection.hpp"
  # include "../oclTraits.hpp"




  /*********************
   ** class oclMatrix **
   ** (declaration)   **
   *********************/
  template <class T>
  class oclMatrix : public Matrix <T>
  {



    public: 


      /**
       * @name            Constructors and destructors
       *                  Constructors and destructors
       */
      //@{

  
      /**
       * @brief           Construct with size 1x1
       */
      inline
      oclMatrix           ()
                        : Matrix <T> (),
                          mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {

        T t;
        Validate (t);
                
      }
      
      
      /**
       * @brief           Construct with dimension array.
       *
       * @param  dim      All 16 dimensions.
       */
      inline
      oclMatrix           (const size_t * dim)
                        : Matrix <T> (dim),
                          mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }


      /**
       * @brief           Construct with dimension and resolution arrays.
       *
       * @param  dim      All 16 dimensions.
       * @param  res      All 16 resolutions.
       */
      inline
      oclMatrix           ( const size_t * dim,
                            const  float * res )
                        : Matrix <T> (dim, res),
                          mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
      
      }
      
      
      /**
       * @brief           Construct 2D (square).
       *
       * @param  n        Dimension of rows and columns.
       */
      inline
      oclMatrix           ( const size_t & n )
                         : Matrix <T> (n),
                           mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }
      
      
      /**
       * @brief           Construct 2D (general).
       *
       * @param  rows     Rows.
       * @param  cols     Columns.
       */
      inline
      oclMatrix           ( const size_t & rows,
                            const size_t & cols )
                         : Matrix <T> (rows, cols),
                           mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
      
      }
      
      
      /**
       * @brief           Construct 3D volume.
       *
       * @param  rows     Rows.
       * @param  cols     Columns.
       * @param  slices   3rd dimension.
       */
      inline
      oclMatrix           ( const size_t & rows,
                            const size_t & cols,
                            const size_t & slices )
                         : Matrix <T> (rows, cols, slices),
                           mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }
      
      
      /**
       * @brief           Construct 4D volume. (or higher dimensional)
       *
       * @param ...       Refer to definition in Matrix <T>.
       */
      inline
      oclMatrix           ( const size_t & col,
                            const size_t & lin,
                            const size_t & cha,
                            const size_t & set,
                            const size_t & eco = 1,
                            const size_t & phs = 1,
                            const size_t & rep = 1,
                            const size_t & seg = 1,
                            const size_t & par = 1,
                            const size_t & slc = 1,
                            const size_t & ida = 1,
                            const size_t & idb = 1,
                            const size_t & idc = 1,
                            const size_t & idd = 1,
                            const size_t & ide = 1,
                            const size_t & ave = 1 )
                         : Matrix <T> (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave),
                           mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }

    
      /**
       * @brief           Transform from Matrix <T> to oclMatrix <T>
       *
       * @param mat       Matrix to copy.
       */
      inline
      oclMatrix           (const Matrix <T> & mat)
                        : Matrix <T> (mat),
                          mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {

        T t;
        Validate (t);
        
      }


      /**
       * @brief           Copy constructor.
       *
       * @param  mat      oclMatrix to copy.
       */
      inline
      oclMatrix           (const oclMatrix <T> & mat)
                        : Matrix <T> ((Matrix<T>) mat),
                          mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size (),
                                                                             * mat.mp_oclData,
                                                                                         true))
      {

        T t;
        Validate (t);

      }

      
      /**
       * @brief           Virtual destructor.
       */
      virtual
      ~oclMatrix          ()
      {

        // delete member oclDataObject (created by oclTraits)
        delete mp_oclData;
    
      }


      //@}
      
      
      /**
       * @name            Elementwise access.
       *                  Elementwise access.
       */
      //@{
      
      
      /**
       * @brief           Copy of p-th element.
       *
       * @param  p        Requested position.
       *
       * @return          Value at p-th position.
       */
      inline
      T
      operator[]          (const size_t & p)
      const;
      
      
      /**
       * @brief           Reference to p-th element.
       *
       * @param  p        Requested position.
       *
       * @return          Value at p-th position.
       */
      inline
      T &
      operator[]          (const size_t & p);
      
      
      //@}
      
      
      /**
       * @brief           assignment operator
       */
      oclMatrix <T> &
      operator=           (const oclMatrix <T> & mat);

      
      /**
       * @name            arithmetic operators
       */
      //@{

      
      /**
       * @brief           Elementwise addition of two matrices.
       *
       * @param mat       Matrix additive.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator+           (const oclMatrix <S> & mat) const;

      
      /**
       * @brief           Elementwise subtraction of two matrices.
       *
       * @param mat       Matrix substruent.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator-           (const oclMatrix <S> & mat) const;
      
      
      /**
       * @brief           Elementwise increment (uniform).
       *
       * @param inc       Increment.
       */
      template <class S>
      oclMatrix <T> &
      operator+=          (const S & inc);


      /**
       * @brief           Elementwise increment (non uniform).
       *
       * @param inc_mat   Matrix containing increments.
       */
      template <class S>
      oclMatrix <T> &
      operator+=          (const oclMatrix <S> & inc_mat);


      /**
       * @brief           Elementwise decrement (uniform).
       *
       * @param dec       Decrement.
       */
      template <class S>
      oclMatrix <T> &
      operator-=          (const S & dec);
      
      
      /**
       * @brief           Elementwise decrement (non uniform).
       *
       * @param dec_mat   Matrix containing decrements.
       */
      template <class S>
      oclMatrix <T> &
      operator-=          (const oclMatrix <S> & dec_mat);

      
      //@}
      
      
      /**
       * @brief           copy relevant data to CPU
       */
      inline
      void
      getData             ()
      const;



    private:


      /**********************
       ** member variables **
       **********************/
             
      // holds data of matrix for calculations on GPU
      oclDataWrapper <T> * mp_oclData;
      
      
      /**********************
       ** static variables **
       **********************/
      
      // verbosity level for operators
      static
      const VerbosityLevel op_v_level = VERB_LOW;
      
      
      /********************************
       ** member functions (private) **
       ********************************/
       
      // allowed element types for an instance of oclMatrix
      void
      Validate            (float  & t)  const {}
      void
      Validate            (size_t & t)  const {}


      
  };




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
      memcpy (Matrix<T> :: _dim, mat.Dim(), INVALID_DIM * sizeof(size_t));
      memcpy (Matrix<T> :: _res, mat.Res(), INVALID_DIM * sizeof( float));

      // temporary save current oclDataObject
      oclDataWrapper <T> * p_old_oclData = this -> mp_oclData;

      /* if buffer is copied, CPU data need not to be copied */
      if (! mat.mp_oclData -> bufferCopyable ())
      {

        // copy data array of Matrix <T>
        Matrix<T> :: _M = mat.Dat();      // uses deep copy of valarray <T>

      }
      
      // create new wrapper for gpu memory (and !new! cpu memory, copy object state!!!)
      /* important to copy object after CPU data (since memory pointer may become invalid) !!! */
      this -> mp_oclData = oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]), Matrix <T> :: Size (), * mat.mp_oclData, true);
      
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
      
    /* copy data to CPU */
    getData ();                     /* TODO: possibility to copy just one element */
      
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

    /* copy data to CPU */
    getData ();                     /* TODO: see above */

    return Matrix <T> :: _M [p];
  
  }



  
  template <class T>
  template <class S>
  inline
  oclMatrix <T>
  oclMatrix <T> ::
  operator+               (const oclMatrix <S> & mat) const
  {

    print_optional ("oclMatrix :: operator+", op_v_level);

    // create matrix for result
    oclMatrix <T> sum (Matrix<T> (this -> Dim()));
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_add (this -> mp_oclData, mat.mp_oclData, sum.mp_oclData, this -> Size ());

    return sum;
    
  }




  template <class T>
  template <class S>
  inline
  oclMatrix <T>
  oclMatrix <T> ::
  operator-               (const oclMatrix <S> & mat) const
  {
    
    print_optional ("oclMatrix :: operator-", op_v_level);
    
    // create matrix for result
    oclMatrix <T> diff (Matrix<T> (this -> Dim()));
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_subtract (this -> mp_oclData, mat.mp_oclData, diff.mp_oclData, this -> Size ());
    
    return diff;
    
  }
  
  
  
  
  template <class T>
  template <class S>
  inline
  oclMatrix <T> &
  oclMatrix <T> ::
  operator-=              (const S & dec)
  {
  
    print_optional ("oclMatrix :: operator-=", op_v_level);
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_dec (this -> mp_oclData, T (dec), this -> Size ());
    
    return *this;
    
  }
  
  
  
  template <class T>
  template <class S>
  inline
  oclMatrix <T> &
  oclMatrix <T> ::
  operator-=              (const oclMatrix <S> & dec_mat)
  {
  
    print_optional ("oclMatrix :: operator-=", op_v_level);
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_subtract (this -> mp_oclData, dec_mat.mp_oclData, this -> mp_oclData, this -> Size ());
    
    return *this;
  
  }
  
  
  

  template <class T>  
  template <class S>
  inline
  oclMatrix <T> &
  oclMatrix <T> ::
  operator+=              (const S & inc)
  {
  
    print_optional ("oclMatrix :: operator+=", op_v_level);
 
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_inc (this -> mp_oclData, T (inc), this -> Size ());
  
    return *this;
  
  }
   
   
   
   
  template <class T>  
  template <class S>
  inline
  oclMatrix <T> &
  oclMatrix <T> ::
  operator+=              (const oclMatrix <S> & inc_mat)
  {
  
    print_optional ("oclMatrix <T> :: operator+=", op_v_level);
 
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_add (this -> mp_oclData, inc_mat.mp_oclData, this -> mp_oclData, this -> Size ());
  
    return *this;
  
  }
      
      
      
  
  template <class T>
  inline
  void
  oclMatrix <T> ::
  getData                 ()
  const
  {
  
    // check if data on CPU are up to date
    if (! this -> mp_oclData -> getSyncState ())
    {
  
      // load GPU data to CPU
      this -> mp_oclData -> getData ();

    }
  
  }




# endif __OCL_MATRIX_HPP__
