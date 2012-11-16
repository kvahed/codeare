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
                                                                   * mat.mp_oclData))
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

      
      //@}
      
      
      /**
       * @brief         copy relevant data to CPU
       */
      inline
      void
      getData           ();



    private:


      /**********************
       ** member variables **
       **********************/
             
      // holds data of matrix for calculations on GPU
      oclDataWrapper <T> * mp_oclData;
      
      
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
      
    std::cout << " ** oclMatrix :: operator=" << std::endl;
      
    // otherwise: self assignment
    if (this != &mat)
    {
        
      // copy members of Matrix <T>
      memcpy (Matrix<T> :: _dim, mat.Dim(), INVALID_DIM * sizeof(size_t));
      memcpy (Matrix<T> :: _res, mat.Res(), INVALID_DIM * sizeof( float));        
      Matrix<T> :: _M = mat.Dat();      // uses deep copy of valarray <T>
      
      // temporary save current oclDataObject
      oclDataWrapper <T> * p_tmp_oclData = mp_oclData;
      
      // create new wrapper for gpu memory (and !new! cpu memory, copy object state!!!)
      mp_oclData = oclTraits <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]), Matrix <T> :: Size (), * mat.mp_oclData);
      
      // delete old oclDataObject
      delete p_tmp_oclData;
      
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
    
    // create matrix for result
    oclMatrix <T> diff (Matrix<T> (this -> Dim()));
    
    // call operator function of oclTraits
    oclTraits <T> :: ocl_operator_subtract (this -> mp_oclData, mat.mp_oclData, diff.mp_oclData, this -> Size ());
    
    return diff;
    
  }
  
  
  

  template <class T>  
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator+=            (const S & inc)
  {
  
    std::cout << "oclMatrix <T> :: operator +=" << std::endl;
  
    // TODO: since calculation is performed on CPU //
    getData ();
  
    // TODO: perform temporarily on CPU //
    for (size_t i = 0; i < this -> Size (); i++)
    {
      this -> _M [i] += T(inc);
    }
    
    // notify modification of CPU data //
    this -> mp_oclData -> setCPUModified ();
  
    return *this;
  
  }
      
      
      
  
  template <class T>
  inline
  void
  oclMatrix <T> ::
  getData               ()
  {
  
    // load GPU data to CPU
    this -> mp_oclData -> getData ();
  
  }




# endif __OCL_MATRIX_HPP__
