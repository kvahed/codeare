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
  # include "../oclDataWrapper.hpp"


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
                      : Matrix <T> ((Matrix<T>) mat)
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
       * @brief         Elementwise addition of two matrices
       *
       * @param mat     Matrix additive.
       */
      template <class S>
      oclMatrix <T>
      operator+         (const oclMatrix <S> & mat) const;
      
      //@}


    private:
      
      oclDataWrapper <T> * mp_oclData;
      
      void
      Validate          (float & t) const {};
      

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
    
    return ((Matrix<T>) *this + (Matrix <S>) mat);
    
  }



# endif __OCL_MATRIX_HPP__
