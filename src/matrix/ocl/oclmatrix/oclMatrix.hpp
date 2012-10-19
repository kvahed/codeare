# ifndef __OCL_MATRIX_HPP__



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
       * @name Constructors and destructors
       *       Constructors and destructors
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
                        mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix<T> :: _M [0]), Matrix<T> :: Size()))
      {
        T t;
        Validate (t);   
        
        oclConnection * ocl = oclConnection :: Instance ();
             
        /* -- */
      }


      /**
       * @brief         Copy constructor
       */
      inline
      oclMatrix         (const oclMatrix <T> & mat)
                      : Matrix <T> ((Matrix<T>) mat),
                        mp_oclData (oclTraits <T> :: make_GPU_Obj (& (Matrix<T> :: _M [0]), Matrix<T> :: Size()))
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
        /* -- */
      }

      //@}
      
      
      /**
       * @name arithmetic operators
       *
       */
      //@{
      
      /**
       * @brief         Elementwise addition of two matrices.
       *
       * @param mat     Matrix additive.
       */
      template <class S>
      oclMatrix <T>
      operator+         (const oclMatrix <S> & mat) const;
      
      /**
       * @brief         Elementwise subtraction of two matrices.
       *
       * @param mat     Matrix substruent.
       */
      template <class S>
      oclMatrix <T>
      operator-         (const oclMatrix <S> & mat) const;
      
      //@}


    private:
      
      // holds data of matrix for calculations on GPU
      oclDataWrapper <T> * mp_oclData;
      
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
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator+             (const oclMatrix <S> & mat) const
  {
    
    oclMatrix <T> sum (Matrix<T> (this -> Dim()));
    
    oclTraits <T> :: ocl_operator_add (this -> mp_oclData, mat.mp_oclData, sum.mp_oclData, this -> Size () );
  
    return sum;
    
  }



  template <class T>
  template <class S>
  oclMatrix <T>
  oclMatrix <T> ::
  operator-             (const oclMatrix <S> & mat) const
  {
    
    oclMatrix <T> diff (Matrix<T> (this -> Dim()));
    
    oclTraits <T> :: ocl_operator_subtract (this -> mp_oclData, mat.mp_oclData, diff.mp_oclData, this -> Size () );
    
    return ((Matrix <T>) *this - (Matrix <T>) mat);
    
  }



# endif __OCL_MATRIX_HPP__
