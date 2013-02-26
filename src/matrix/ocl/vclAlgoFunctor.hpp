# ifndef __VCL_ALGO_FUNCTOR_HPP__




  /************
   ** makros **
   ************/
  # define __VCL_ALGO_FUNCTOR_HPP__



  /**************
   ** includes **
   **************/
   
  // ocl
  # include "oclSettings.hpp"
  
  // ViennaCl
  # include "/usr/include/viennacl/linalg/prod.hpp"




  /**************************
   ** forward declarations **
   **************************/
  template <class T, class S>
  class oclViennaClObject;




  /***************************
   ** class: vclAlgoFunctor **
   **   (interface)         **
   **   (still containing   **
   **    dummy functions)   **
   ***************************/
  template <class T, class S>
  class vclAlgoFunctor
  {
  
  

    public:


      /**
       * @brief               overloaded operator (functor)
       */
      inline
      virtual
      void
      operator()              ()
      const
      {
      
        std::cout << "not implemented" << std::endl;
        
      }
    
    

    protected:
    
    
      /**
       * @name                constructors and destructors
       */
      //@{
      
      
      /**
       * @brief               constructor
       */
      vclAlgoFunctor          ( const oclViennaClObject <T, S> * const p_vclObj )
                             : mp_vclObj (p_vclObj)
      {
      
        print_optional ("Ctor: \"vclAlgoFunctor\"", VERB_HIGH);
        
      }
      
      
      //@}
      
      
      /**********************
       ** member variables **
       **********************/
      
      const oclViennaClObject <T, S> * const mp_vclObj;


      /**********************
       ** static variables **
       **********************/
      static
      const VerbosityLevel op_v_level = VERB_MIDDLE;



  }; // class vclAlgoFunctor
  
  

  
  /****************************
   ** class: vclSubtractAlgo **
   **    (derived)           **
   ****************************/
  template <class T, class S>
  class vclAddSubAlgo : public vclAlgoFunctor <T, S>
  {

  
  
    public:

    
      /**
       * @name              constructors and destructors
       */
      //@{
      
      
      /**
       * @brief             constructor
       */
      vclAddSubAlgo       ( const oclViennaClObject <T, S> * const p_vclObj,
                            const               int                 sign )
                           : vclAlgoFunctor <T, S> (p_vclObj),
                             m_sign                    (sign)
      {

        print_optional ("Ctor: \"vclAddSubAlgo\"", VERB_HIGH);

      }
      
      
      //@}
      
      
      /**
       * @brief             refer to base class
       */
      inline
      virtual
      void
      operator()            ()
      const
      {

        print_optional ("vclAddSubAlgo <T, S> :: operator()", vclAlgoFunctor <T, S> :: op_v_level);

        /**
         * create ViennaCl arguments
         */
        viennacl :: vector <T>   arg1 = vclAlgoFunctor <T, S> :: mp_vclObj -> template getVCLVector <T> (0);
        viennacl :: vector <S>   arg2 = vclAlgoFunctor <T, S> :: mp_vclObj -> template getVCLVector <S> (1);
        viennacl :: vector <T> result = vclAlgoFunctor <T, S> :: mp_vclObj -> template getVCLVector <T> (2);
        
        /* perform calculation */
        if (m_sign < 0)
          result = arg1 - arg2;
        else
          result = arg1 + arg2;
        
      }
      
      
    private:
    
      int m_sign;
  
  
  
  }; // class vclSubtractAlgo



  /****************************
   ** class: vclMatProdAlgo  **
   **    (derived)           **
   ****************************/
  template <class T, class S>
  class vclMatProdAlgo : public vclAlgoFunctor <T, S>
  {

  
  
    public:

    
      /**
       * @name              constructors and destructors
       */
      //@{
      
      
      /**
       * @brief             constructor
       */
      vclMatProdAlgo      ( const oclViennaClObject <T, S> * const p_vclObj )
                           : vclAlgoFunctor <T, S> (p_vclObj),
                                          m     (p_vclObj -> getScalarArg (0)),
                                          k     (p_vclObj -> getScalarArg (1)),
                                          n     (p_vclObj -> getScalarArg (2)),
                                     transA     (p_vclObj -> getScalarArg (3)),
                                     transB     (p_vclObj -> getScalarArg (4))
      {

        print_optional ("Ctor: \"vclMatProdAlgo\"", VERB_HIGH);

      }
      
      
      //@}
      
      
      /**
       * @brief             refer to base class
       */
      inline
      virtual
      void
      operator()            ()
      const
      {

        print_optional ("vclMatProdAlgo <T> :: operator()", vclAlgoFunctor <T, S> :: op_v_level);

        /**
         * create ViennaCl arguments
         */
        viennacl :: matrix <T, viennacl :: column_major>   arg1 = vclAlgoFunctor <T, S> :: mp_vclObj
                                                                   -> template getVCLMatrix <T, viennacl :: column_major> (0, m, k);
        viennacl :: matrix <T, viennacl :: column_major>   arg2 = vclAlgoFunctor <T, S> :: mp_vclObj
                                                                   -> template getVCLMatrix <S, viennacl :: column_major> (1, k, n);
        viennacl :: matrix <T, viennacl :: column_major> result = vclAlgoFunctor <T, S> :: mp_vclObj
                                                                   -> template getVCLMatrix <T, viennacl :: column_major> (2, m, n);
        
        using viennacl :: trans;
        
        /* perform calculation */
        if (transA)
          if (transB) result = viennacl :: linalg :: prod (trans (arg1), trans (arg2));
          else        result = viennacl :: linalg :: prod (trans (arg1),        arg2 );
        else
          if (transB) result = viennacl :: linalg :: prod (       arg1 , trans (arg2));
          else        result = viennacl :: linalg :: prod (       arg1 ,        arg2 );
        
        viennacl::ocl::get_queue ().finish ();

      }
    
    
    private:
      
      int m, k, n;
      int transA, transB;
  
  
  }; // class vclMatProdAlgo

  
  
  
  /**********************
   ** global functions **
   **********************/
  template <class T, class S>
  static
  const vclAlgoFunctor <T, S> * const
  get_algo_functor     (const       vclAlgoType     &           algo,
                        const oclViennaClObject <T, S> * const p_vclObj)
  {
        
    print_optional (" :: get_algo_functor", VERB_HIGH);

    /* choose requested algorithm */
    switch (algo)
    {
    
      case vclSUBTRACT:
        return new vclAddSubAlgo <T, S> (p_vclObj, -1);
      
      case vclMATPROD:
        return new vclMatProdAlgo <T, S> (p_vclObj);

      default:
        throw "*!* Specified ViennaCL algorithm not available! *!*";
    
    }

  }
  
  
  
  
# endif // __VCL_ALGO_FUNCTOR_HPP__
