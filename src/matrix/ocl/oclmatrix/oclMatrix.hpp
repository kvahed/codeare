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
       * @name          Constructors and destructors
       *                Constructors and destructors
       */
      //@{

  
      /**
       * @brief         Construct with size 1x1
       */
      inline
      oclMatrix         ()
                      : Matrix <T> (),
                        mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix<T> :: _M [0]), Matrix<T> :: Size()))
      {

        T t;
        Validate (t);
                
        /* -- */
      };

    
      /**
       * @brief         Transform from Matrix <T> to oclMatrix <T>
       */
      inline
      oclMatrix         (const Matrix <T> & mat)
                      : Matrix <T> (mat),
                        mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                   Matrix <T> :: Size ()))
      {

        T t;
        Validate (t);   
                     
        /* -- */
      }


      /**
       * @brief         Copy constructor
       */
      inline
      oclMatrix         (const oclMatrix <T> & mat)
                      : Matrix <T> ((Matrix<T>) mat),
                        mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                   Matrix <T> :: Size (),
                                                                   * mat.mp_oclData,
                                                                   true))
      {

        T t;
        Validate (t);

        /* -- */
      }

      
      /**
       * @brief         Virtual destructor
       */
      virtual
      ~oclMatrix        ()
      {

        // delete member oclDataObject (created by oclTraits)
        delete mp_oclData;
    
      }


      //@}
      
      
      /**
       * @brief         assignment operator
       */
      oclMatrix <T> &
      operator=         (const oclMatrix <T> & mat);

      
      /**
       * @name          arithmetic operators
       */
      //@{

      
      /**
       * @brief         Elementwise addition of two matrices.
       *
       * @param mat     Matrix additive.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator+         (const oclMatrix <S> & mat) const;

      
      /**
       * @brief         Elementwise subtraction of two matrices.
       *
       * @param mat     Matrix substruent.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator-         (const oclMatrix <S> & mat) const;
      
      
      /**
       * @brief         Elementwise increment.
       *
       * @param inc     Increment.
       */
      template <class S>
      oclMatrix <T>
      operator+=        (const S & inc);


      /**
       * @brief         Elementwise increment.
       *
       * @param inc     Increment.
       */
      template <class S>
      oclMatrix <T>
      operator+=        (const oclMatrix <S> & inc_mat);


      /**
       * @brief         Elementwise decrement.
       *
       * @param dec     Decrement.
       */
      template <class S>
      oclMatrix <T>
      operator-=        (const S & dec);

      
      //@}
      
      
      /**
       * @brief         copy relevant data to CPU
       */
      inline
      void
      getData           ()
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
      const VerbosityLevel op_v_level = VERB_MIDDLE;
      
      
      /********************************
       ** member functions (private) **
       ********************************/
       
      // allowed element types for an instance of oclMatrix
      void
      Validate          (float  & t)  const {}
      void
      Validate          (size_t & t)  const {}


      
  };




  /******************************
   ** definitions in oclMatrix **
   ******************************/



  template <class T>
  oclMatrix <T> &
  oclMatrix <T> ::
  operator=            (const oclMatrix <T> & mat)
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



  
  template <class T>
  template <class S>
  inline
  oclMatrix <T>
  oclMatrix <T> ::
  operator+             (const oclMatrix <S> & mat) const
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
  operator-             (const oclMatrix <S> & mat) const
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
  oclMatrix <T>
  oclMatrix <T> ::
  operator-=            (const S & dec)
  {
  
    print_optional ("oclMatrix :: operator-=", op_v_level);
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_dec (this -> mp_oclData, T (dec), this -> Size ());
    
    return *this;
    
  }
  
  
  

  template <class T>  
  template <class S>
  inline
  oclMatrix <T>
  oclMatrix <T> ::
  operator+=            (const S & inc)
  {
  
    print_optional ("oclMatrix :: operator+=", op_v_level);
 
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_inc (this -> mp_oclData, T (inc), this -> Size ());
  
    return *this;
  
  }
   
   
   
   
  template <class T>  
  template <class S>
  inline
  oclMatrix <T>
  oclMatrix <T> ::
  operator+=            (const oclMatrix <S> & inc_mat)
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
  getData               ()
  const
  {
  
    // load GPU data to CPU
    this -> mp_oclData -> getData ();
  
  }




# endif __OCL_MATRIX_HPP__
